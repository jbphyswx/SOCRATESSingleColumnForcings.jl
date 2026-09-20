#=
Regridding pipeline: build/evaluate per-column interpolants along a dimension (`interp_along_dim`,
`var_to_new_coord`), and the unified `regrid_to_z_and_time` that regrids a source onto a new vertical
grid and then builds time splines. Source-specific steps are generic verbs dispatched on
`AbstractRegridSource` (no `isa` branching in the shared body).

The interpolation API lives in the `Interpolation` submodule; call it via qualified `Interpolation.foo`.
=#

"""
    ConservativeRegridingOpts(; k = 1, f_enhancement_factor = 1, f_p_enhancement_factor = 1,
                     kwargs = Interpolation.default_conservative_interp_kwargs)

Settings for a conservative vertical regrid: the refinement factors the mass matrix is built at, and
the solver options in [`Interpolation.default_conservative_interp_kwargs`](@ref).

Its presence is what selects a conservative regrid — a stage holds one of these or `nothing`, so
there is no separate flag that can disagree with the settings.
"""
struct ConservativeRegridingOpts{CK <: Interpolation.DCIKT}
    k::Int
    f_enhancement_factor::Int
    f_p_enhancement_factor::Int
    kwargs::CK
end
ConservativeRegridingOpts(;
    k::Int = 1,
    f_enhancement_factor::Int = 1,
    f_p_enhancement_factor::Int = 1,
    kwargs::Interpolation.DCIKT = Interpolation.default_conservative_interp_kwargs,
) = ConservativeRegridingOpts(k, f_enhancement_factor, f_p_enhancement_factor, kwargs)

"""
    RegriddingOpts(; method, bc, conservative, drop_collinear)

How a field is evaluated onto the new vertical grid. `conservative` is a [`ConservativeRegridingOpts`](@ref)
or `nothing`.

`bc` defaults to extrapolation when conservative, because the regrid integrates past the outermost
knots.
"""
struct RegriddingOpts{
    M <: Interpolation.AbstractInterpolationMethod,
    BC <: Interpolation.ValidBoundaryConditions,
    C <: Union{Nothing, ConservativeRegridingOpts},
    DC <: Val,
}
    method::M
    bc::BC
    conservative::C
    drop_collinear::DC
end
RegriddingOpts(;
    method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
    conservative::Union{Nothing, ConservativeRegridingOpts} = nothing,
    bc::Interpolation.ValidBoundaryConditions = isnothing(conservative) ? Interpolation.ErrorBoundaryCondition() : Interpolation.ExtrapolateBoundaryCondition(),
    drop_collinear::Val = Val(false),
) = RegriddingOpts(method, bc, conservative, drop_collinear)


"""
Cache key for a conservative mass matrix: everything the matrix is a function of — the interpolation
method, the boundary condition, the refinement factor, and the target grid.
"""
const MassMatrixKey = Tuple{DataType, DataType, Int, <:Vector{<:AbstractFloat}}

"""
    MassMatrixCache{FT}()

Conservative mass matrices and their factorizations, reusable across fields and across calls.

Keyed by [`MassMatrixKey`](@ref)
"""
struct MassMatrixCache{FT}
    A::Dict{MassMatrixKey, Matrix{FT}}
    Af::Dict{MassMatrixKey, LinearAlgebra.LU{FT, Matrix{FT}, Vector{Int}}}
end
MassMatrixCache{FT}() where {FT} = MassMatrixCache{FT}(
    Dict{MassMatrixKey, Matrix{FT}}(),
    Dict{MassMatrixKey, LinearAlgebra.LU{FT, Matrix{FT}, Vector{Int}}}(),
)

"""
    mass_matrix!(cache, z_new, opts) -> (A, Af)

The conservative mass matrix for `z_new` under `opts` and its factorization, building and storing them
on first use. Returns `(nothing, nothing)` when `opts` is not conservative.
"""
function mass_matrix!(cache::MassMatrixCache{FT}, z_new, opts::RegriddingOpts) where {FT}
    isnothing(opts.conservative) && return (nothing, nothing)
    key = (typeof(opts.method), typeof(opts.bc), opts.conservative.k, FT.(collect(z_new)))::MassMatrixKey
    if !haskey(cache.A, key)
        cache.A[key] = Matrix{FT}(Interpolation.conservative_mass_matrix(
            z_new; method = opts.method, k = opts.conservative.k, bc = opts.bc,
        ))
        cache.Af[key] = LinearAlgebra.factorize(cache.A[key])
    end
    return (cache.A[key], cache.Af[key])
end

"""
    inplace_eval_kernel(method, conservative)

The in-place `(y, x, xp, fp; bc)` evaluation kernel for `method`, or `nothing` when it has none and the
caller must use the allocating form. A conservative regrid solves a system per column and has no
in-place form, so it returns `nothing` for every method.
"""
inplace_eval_kernel(::Interpolation.AbstractInterpolationMethod, ::ConservativeRegridingOpts) = nothing
inplace_eval_kernel(::Interpolation.AbstractInterpolationMethod, ::Nothing) = nothing
inplace_eval_kernel(::Interpolation.FastLinear1DInterpolationMethod, ::Nothing) = Interpolation.fast1d_interp!

"""
    _alloc_eval_output(vardata, interp_dim_num, x_new, xp)

Output array for evaluating `vardata`'s columns at `x_new`: `vardata`'s shape with `interp_dim_num`
resized to `length(x_new)`
"""
_alloc_eval_output(vardata, interp_dim_num::Int, x_new, xp) = Array{
    promote_type(eltype(x_new), eltype(xp), eltype(vardata)),
}(
    undef,
    ntuple(d -> d == interp_dim_num ? length(x_new) : size(vardata, d), Val(ndims(vardata)))...,
)

"""
    eval_columns_into!(out, interp_dim_num, x_new, columns, interp_func, kernel, bc)

Evaluate each `(coordinate, values)` pair of `columns` at `x_new` into `out`, in the order `eachslice`
visits it.
"""
function eval_columns_into!(out, interp_dim_num::Int, x_new, columns, interp_func, kernel, bc)
    column_dims = Tuple(d for d in 1:ndims(out) if d != interp_dim_num)
    for (dst, (xp, fp)) in zip(eachslice(out; dims = column_dims), columns)
        if isnothing(kernel)
            dst .= interp_func(x_new, xp, fp)
        else
            kernel(dst, x_new, xp, fp; bc = bc)
        end
    end
    return out
end

"""
    interp_along_dim(var, interp_dim, interp_dim_in,
                     interpolant_coord_types, interpolant_value_types, drop_collinear_val;
                     interp_dim_out, bc, method, conservative, ...)

Interpolate `var` along `interp_dim` (a dimension number, or a name resolved against `data`), from
the coordinate `interp_dim_in` it currently sits on.

`interp_dim_out` selects between the two paths:
- given, it EVALUATES at those coordinates and returns an array;
- `nothing`, it BUILDS one interpolant per column and returns them, storing nodes and values in the
  backings named by `interpolant_coord_types` / `interpolant_value_types` and pruning collinear
  interior nodes per `drop_collinear_val`.

`interp_dim_in_is_full_array` says whether `interp_dim_in` is a coordinate per column (same shape as
`var`) rather than one shared template vector.

`bc` is the out-of-range policy

`conservative` is a [`ConservativeRegridingOpts`](@ref) or `nothing`;
"""
function interp_along_dim(
    var,
    interp_dim::Union{String, Int},
    interp_dim_in,
    ::Type{interpolant_coord_types} = Tuple{Vector, Nothing},   # (backing, eltype) spec TYPE for the coordinate; `Nothing` = passthrough
    ::Type{interpolant_value_types} = Tuple{Vector, Nothing},   # (backing, eltype) spec TYPE for the values
    drop_collinear_val::Val{drop_collinear} = Val(false),
    ;
    interp_dim_out = nothing,
    data = nothing,
    data_func::Union{Function, Nothing} = nothing,
    interp_dim_in_is_full_array::Bool = true,
    # reshape_ground::Bool = true,
    method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
    conservative::Union{Nothing, ConservativeRegridingOpts} = nothing,
    # A conservative regrid has to integrate past the outermost knots, so it extrapolates by default.
    bc::Interpolation.ValidBoundaryConditions = isnothing(conservative) ? Interpolation.ErrorBoundaryCondition() : Interpolation.ExtrapolateBoundaryCondition(),
    A::Union{Nothing, AbstractArray} = nothing,
    Af::Union{Nothing, AbstractArray, LinearAlgebra.Factorization} = nothing, # precomputed factorization of A :: technically, A being Diagonal, or Triangular or something could lead to AbstractMatrix Af so we allow both AbstractMatrix and Factorization
) where {interpolant_coord_types <: Tuple, interpolant_value_types <: Tuple, drop_collinear}



    # if data is a string, read the data out from data (creates data `vardata` from `var` whether var is string or already is data)
    vardata = isa(var, String) ? data[var] : var


    # get the interpolation dimension and combine air and ground data
    interp_dim_num = dim_num(interp_dim, vardata)
    # Materialize any lazy DiskArray (e.g. BroadcastDiskArray from NCDatasets arithmetic) and
    # strip Union{Missing,Float64} → Float64.  Must happen AFTER dim_num while labels exist.
    vardata = _materialize(vardata)
    # vardata = combine_air_and_ground_data(vardata, vardatag, interp_dim_num; ...) # (unused)

    if !isnothing(data_func) # apply data_func if we need to
        vardata = data_func(vardata)
    end

    # Unified evaluator with a positional (xnew, xp, fp) interface: conservative regridder or plain interpolation.
    interp_func =
        isnothing(conservative) ?
        ((xnew, xp, fp) -> Interpolation.interpolate_1d(xnew, xp, fp, method; bc = bc)) :
        ((xnew, xp, fp) -> Interpolation.conservative_regridder(
            xnew,
            xp,
            fp;
            method = method,
            bc = bc,
            k = conservative.k,
            f_enhancement_factor = conservative.f_enhancement_factor,
            f_p_enhancement_factor = conservative.f_p_enhancement_factor,
            A = A,
            Af = Af,
            conservative.kwargs...,
        ))
    # coordinate and value get INDEPENDENT storage specs (backing, eltype); Base.Fix1 fixes the spec TYPE, varies v.
    coord_conv = Base.Fix1(Interpolation.coerce_vector, interpolant_coord_types)
    value_conv = Base.Fix1(Interpolation.coerce_vector, interpolant_value_types)

    # mapslices to apply along timedim, see https://docs.julialang.org/en/v1/base/arrays/#Base.mapslices
    if !interp_dim_in_is_full_array
        if isnothing(interp_dim_out)
            # BUILD path: one interpolant per column, pruned per `drop_collinear` during construction
            # (no full-node interpolant materialized). `map` takes the array eltype from the results —
            # concrete when the backing type is length-invariant (Vector, …), abstract when it encodes

            # One path for any AbstractVector backing and either drop setting.
            xin = coord_conv(interp_dim_in)
            result = map(eachslice(vardata; dims = Tuple(d for d in 1:ndims(vardata) if d != interp_dim_num))) do d
                Interpolation.build_spline(method, xin, value_conv(d); bc = bc, drop_collinear = drop_collinear_val)
            end
            return add_dim(result, interp_dim_num) # reinsert the size-1 interp dim to match the evaluate-path layout

        else
            # Every column shares this coordinate, so coerce it once rather than once per column.
            xin_shared = coord_conv(interp_dim_in)
            cols = eachslice(vardata; dims = Tuple(d for d in 1:ndims(vardata) if d != interp_dim_num))
            out = _alloc_eval_output(vardata, interp_dim_num, interp_dim_out, xin_shared)
            return eval_columns_into!(
                out, interp_dim_num, interp_dim_out,
                ((xin_shared, value_conv(d)) for d in cols),
                interp_func, inplace_eval_kernel(method, conservative), bc,
            )
        end
    else # vectorize over input dim values as well as data (no support for vectorize over output dim yet)
        if isnothing(interp_dim_out)
            # BUILD path (stacked [coord data] layout): same unified map as the non-full-array case;
            # each slice is the 2-column matrix (d[:,1] = coord, d[:,2] = data), pruned per
            # `drop_collinear` during construction; `map` sets the eltype from the results. Only this
            # path needs the stacked array; the evaluate path below walks the two inputs together.
            catd = ndims(vardata) + 1
            _input = cat(interp_dim_in, vardata; dims = catd) # although maybe the input is just a vector in which case this won't work..... we could just pass it in
            result = map(eachslice(_input; dims = Tuple(d for d in 1:ndims(_input) if d != interp_dim_num && d != catd))) do d
                Interpolation.build_spline(method, coord_conv(d[:, 1]), value_conv(d[:, 2]); bc = bc, drop_collinear = drop_collinear_val)
            end
            return add_dim(result, interp_dim_num) # reinsert the size-1 interp dim (catd reduced out)
        else

            column_dims = Tuple(d for d in 1:ndims(vardata) if d != interp_dim_num)
            out = _alloc_eval_output(vardata, interp_dim_num, interp_dim_out, interp_dim_in)
            return eval_columns_into!(
                out, interp_dim_num, interp_dim_out,
                (
                    (coord_conv(xc), value_conv(fc)) for (xc, fc) in
                    zip(eachslice(interp_dim_in; dims = column_dims), eachslice(vardata; dims = column_dims))
                ),
                interp_func, inplace_eval_kernel(method, conservative), bc,
            )
        end
    end
end

"""
    var_to_new_coord(var, coord_in, interp_dim, interpolant_coord_types, interpolant_value_types,
                     drop_collinear; coord_new, weight, weight_regridded, bc, ...)

[`interp_along_dim`](@ref) with optional weighting: `coord_new = nothing` returns built interpolants,
otherwise it evaluates onto `coord_new`.

`weight` makes the regrid mass-weighted, computing `interp(var .* weight) ./ interp(weight)` — for an
extensive variable, or for conservative regridding, this is typically density. The denominator is the
same for every field sharing a `(weight, grid, method)`, so a caller that regrids many fields can
compute it once and pass it as `weight_regridded` instead of having it recomputed per field.
"""
function var_to_new_coord(
    var,
    coord_in,
    interp_dim::Union{String, Number},
    ::Type{interpolant_coord_types} = Tuple{Vector, Nothing},   # (backing, eltype) spec TYPE for the coordinate
    ::Type{interpolant_value_types} = Tuple{Vector, Nothing},   # (backing, eltype) spec TYPE for the values
    ::Val{drop_collinear} = Val(false), # prune collinear nodes of built interpolants (build path only)
    ;
    coord_new = nothing,
    data = nothing,
    data_func::Union{Function, Nothing} = nothing,
    method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
    conservative::Union{Nothing, ConservativeRegridingOpts} = nothing,
    bc::Interpolation.ValidBoundaryConditions = isnothing(conservative) ? Interpolation.ErrorBoundaryCondition() : Interpolation.ErrorBoundaryCondition(),
    weight = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    # Precomputed regridded weight = the `interp(weight)` normalization denominator. When supplied it is
    # used directly, skipping the redundant second `interp_along_dim(weight, …)`: that denominator is
    # identical for every field regridded with the same (weight, grid, method), so the caller computes it once.
    weight_regridded = nothing,
    A::Union{Nothing, AbstractArray} = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    Af::Union{Nothing, AbstractArray, LinearAlgebra.Factorization} = nothing, # precomputed factorization of A :: technically, A being Diagonal, or Triangular or something could lead to AbstractMatrix Af so we allow both AbstractMatrix and Factorization
) where {interpolant_coord_types <: Tuple, interpolant_value_types <: Tuple, drop_collinear}
    vardata = isa(var, String) ? data[var] : var
    if ~isnothing(data)
        interp_dim = isa(interp_dim, String) ? dim_num(interp_dim, data) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
    else
        interp_dim = isa(interp_dim, String) ? dim_num(interp_dim, vardata) : interp_dim # if interp_dim is a string, you need to provide the underlying data so we can get this dimension
    end

    # Materialize any lazy DiskArray and strip Union{Missing,Float64} before vardata .* weight
    # which would otherwise create a new BroadcastDiskArray.  Must happen AFTER interp_dim
    # resolution while dimension labels may still be present on vardata.
    vardata = _materialize(vardata)
    isnothing(weight) || (weight = _materialize(weight))

    # evaluate interp_dim_in_is_full_array based on the size of the input... interp_dim_in is full array false is much faster cause dont have to double loop in vectorization...

    interp_dim_in_is_full_array = (size(coord_in) == size(vardata))
    _interp(v) = interp_along_dim(
        v,
        interp_dim,
        coord_in,
        interpolant_coord_types,
        interpolant_value_types,
        Val(drop_collinear);
        interp_dim_out = coord_new,
        data = data,
        data_func = data_func,
        interp_dim_in_is_full_array = interp_dim_in_is_full_array,
        method = method,
        conservative = conservative,
        bc = bc,
        A = A,
        Af = Af,
    )

    isnothing(weight) && return _interp(vardata)
    # Denominator = interp(weight). Reuse the caller-precomputed value if given (identical across every
    # field sharing this weight/grid/method); otherwise compute it here.
    den = isnothing(weight_regridded) ? _interp(weight) : weight_regridded
    return _interp(vardata .* weight) ./ den
end


# ================================================================================================ #
# Regrid sources
#
# `regrid_to_z_and_time` runs ONE pipeline (regrid to `z_new`, then build time splines) across data
# sources. The handful of source-specific steps are generic verbs dispatched on the source type
# (no `isa` branching in the shared body). Add a source by subtyping `AbstractRegridSource` and
# defining the verbs below for it.
# ================================================================================================ #
abstract type AbstractRegridSource end

"ERA5-derived Atlas forcing *input* files (the default source)."
struct AtlasInput <: AbstractRegridSource end

"CliMA LES *output* files: already on a z grid; time stored in days from the run start."
struct LESOutput{FT <: AbstractForcingType} <: AbstractRegridSource
    forcing_type::FT
end

# dataset to read from when only a variable name is given (Atlas always passes `data`; LES self-opens)
regrid_source_data(::AbstractRegridSource, flight_number, data) = data
regrid_source_data(source::LESOutput, flight_number, data) =
    isnothing(data) ? open_atlas_les_output(flight_number, source.forcing_type).data : data

# old vertical coordinate (before regridding onto `z_new`)
regrid_source_z_old(::AtlasInput, data; thermodynamics_backend) =
    z_from_data(data; thermodynamics_backend) # time-varying; ground value pads the bottom with 0
regrid_source_z_old(::LESOutput, data, flight_number::Int, ::Type{FT}; thermodynamics_backend) where {FT} =
    # vec(Array(data["z"])) # static column, already present in the LES output
    _load_grid(flight_number, Val(true), FT) # static column # is float64 so less lossy

# old time coordinate, in SECONDS from the start of the data
regrid_source_t_old(::AtlasInput, data) = vec(Array(data["tsec"]))
"""
Because this is stored at float32 precision but is days, it has some jitter.
It appears to be data every 300s. data loooks like:

Float32[
 36.001736
 36.005207
 36.008682
 36.012154
 ...
 36.491318
 36.494793
 36.498264
 ]

 You can verify with:

    t = vec(Array(SSCF.open_atlas_les_output(9, SSCF.ObsForcing()).data["time"]))
    n = length(t)
    Δ = 300/86400
    ideal = Float32[Float32(Float64(t[1]) + k*Δ) for k in 0:n-1]   # true grid point, rounded to Float32 ONCE
    all(prevfloat.(ideal) .<= t .<= nextfloat.(ideal))             # => true

that it's all rounding error...
"""
function regrid_source_t_old(::LESOutput, data, ::Val{fix_rounding_error} = Val(true), ::Val{fix_averaging} = Val(true)) where {fix_rounding_error, fix_averaging} # Val for type stability
    if fix_rounding_error
        n = length(data["time"])
        out = fix_averaging ? Vector{Int64}(undef, n+1) : Vector{Int64}(undef, n) # could do like UInt32 but not worth extra readability loss
        t = vec(Array(data["time"]))
        @inbounds for i in eachindex(fix_averaging ? view(out, 1:n) : out)
            out[i] = Int64(300 * round((t[i] - t[1]) * 86400 / 300))
            if !fix_averaging # data remains on midpoints
                out[i] += Int64(150) # works on sliced data, etc
            end
        end
        if fix_averaging
            out[end] = out[end-1] + Int64(300) # add the last point for the final window's end
        end
    else
        t = vec(Array(data["time"])) # LES stores time in days
        # out = @. (t - t[1]) * (24 * 3600) # -> seconds from start
        # out = @. (t - t[1]) * (24 * 3600) # -> seconds from start
        if fix_averaging
            out = Vector{eltype(t)}(undef, length(t)+1)
            @. out[1:end-1] = (round((t - t[1]) * (24 * 3600)))
            out[end] = out[end-1] + Int64(300) # 
        else
            out = @. (t - t[1]) * (24 * 3600) + Int64(150) # seconds from start, data remains on midpoints
        end
    end
    return out
end

"""
    les_output_period(t_days)

Sampling period [s] of an LES time axis given in days
"""
@inline les_output_period(::Type{T} = Int64) where {T} =  T(300)

"""
    les_averaging_window(t_days) -> (; period, half_period)

The averaging window each LES sample represents, in seconds: its length, and the offset from a
sample's label back to the start of its window.

Every value in an LES output file is a mean over one window, labelled at that window's MIDPOINT. SAM
accumulates over `nstat` steps between writes and stamps the record `day - nstat*dt/2/86400`
"""
@inline les_averaging_window(::Type{T} = Int64) where {T} =  (; period = T(300), half_period = T(150))

# index of the reference/initial timestep within `t_old`
regrid_source_initial_ind(::AtlasInput, data, flight_number, t_old) =
    initial_index(data, flight_number; t_old = t_old)
regrid_source_initial_ind(::LESOutput, data, flight_number, t_old) = 1 # LES output starts at the initial condition

# Instant that model time 0 refers to
regrid_source_t_origin(::AbstractRegridSource, data, t_old, initial_ind) = t_old[initial_ind] # `tsec` is absolute (seconds after 00Z on nbdate)
regrid_source_t_origin(::LESOutput, data, t_old, initial_ind) = zero(eltype(t_old)) # `regrid_source_t_old` already opens the axis at the first window's midpoint

# is `z_old` time-varying (needs per-timestep selection) rather than a single static column?
regrid_source_z_time_varying(::AtlasInput) = true
regrid_source_z_time_varying(::LESOutput) = false

# must the z axis be reversed to run ground -> top?
regrid_source_reverse_z(::AtlasInput) = true  # Atlas pressure levels come top -> ground
regrid_source_reverse_z(::LESOutput) = false  # LES output is already ground -> top

"""
    regrid_to_z_and_time(var, z_new, z_dim, time_dim, flight_number,
                         interpolant_coord_types, interpolant_value_types;
                         source = AtlasInput(), bc, kwargs...)

Take data from a regrid `source`, interpolate it onto `z_new`, then build time splines over the times
we have. To vectorize properly over `z_new` it should be the same shape as `vardata` (+ ground). The
`source` (default [`AtlasInput`](@ref); [`LESOutput`](@ref) for LES files) selects the few
source-specific steps via the `regrid_source_*` verbs; everything else is shared.

The two stages are configured independently: `z_regrid_opts` is a [`RegriddingOpts`](@ref) for the
evaluation onto `z_new`, `t_regrid_opts` a [`RegriddingOpts`](@ref) for the returned time splines.
`interpolant_coord_types` / `interpolant_value_types` are the `Tuple{Backing, Eltype}` storage specs
for the built interpolants' nodes and values.
"""
function regrid_to_z_and_time(
    var,
    z_new,
    z_dim::Union{String, Number},
    time_dim::Union{String, Number},
    flight_number::Int,
    ::Type{interpolant_coord_types} = Tuple{StepRangeLen, Nothing},  # stored time-axis coordinate: a range (O(1) eval), eltype as-read
    ::Type{interpolant_value_types} = Tuple{Vector, Float64},        # values (and the transient z coordinate): backing + eltype
    ;
    source::AbstractRegridSource = AtlasInput(),
    thermodynamics_backend = DefaultThermodynamicsBackend(),
    varg = nothing,
    z_old = nothing,
    t_old = nothing,
    data = nothing,
    initial_condition::Bool = false,
    assume_monotonic::Bool = false,
    z_regrid_opts::RegriddingOpts = RegriddingOpts(),
    t_regrid_opts::RegriddingOpts = RegriddingOpts(),
    ground_indices = :end,
    weight = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    weightg = nothing, # for extensive variables and for conservative regridding, you may wish to weight by something like density when interpolating in z...
    return_before_interp::Bool = false,
    A::Union{Nothing, AbstractArray} = nothing,
    Af::Union{Nothing, AbstractArray, LinearAlgebra.Factorization} = nothing, # precomputed factorization of A :: technically, A being Diagonal, or Triangular or something could lead to AbstractMatrix Af so we allow both AbstractMatrix and Factorization
    # Precomputed regridded-weight denominator `interp(weight)` (see `var_to_new_coord`), passed straight
    # through to the z-interpolation so the caller can compute it ONCE and reuse it across fields.
    weight_regridded = nothing,
    # Return the z-interpolated field immediately, before building time splines. The caller uses this to
    # obtain the regridded-weight denominator once (an unweighted regrid of ρ) for reuse across fields.
    return_after_z_interp::Bool = false,
) where {interpolant_coord_types <: Tuple, interpolant_value_types <: Tuple}

    if !isnothing(z_regrid_opts.conservative) && isnothing(A)
        A = Interpolation.conservative_mass_matrix(
            z_new;
            method = z_regrid_opts.method,
            bc = z_regrid_opts.bc,
            k = z_regrid_opts.conservative.k,
        )
        Af = LinearAlgebra.factorize(A)
    end

    # Preprocess the weights by re-running with the weights as the data, short-circuiting via
    # `return_before_interp` (`varg = weightg` covers both the ground and no-ground cases).
    if !isnothing(weight)
        weight = regrid_to_z_and_time(
            weight,
            z_new,
            z_dim,
            time_dim,
            flight_number,
            interpolant_coord_types,
            interpolant_value_types;
            source = source,
            weight = nothing,
            weightg = nothing,
            varg = weightg,
            return_before_interp = true,
            thermodynamics_backend = thermodynamics_backend,
            z_old = z_old,
            t_old = t_old,
            data = data,
            initial_condition = initial_condition,
            assume_monotonic = assume_monotonic,
            z_regrid_opts = z_regrid_opts,
            t_regrid_opts = t_regrid_opts,
            ground_indices = ground_indices,
            A = A,
            Af = Af,
        )
    end

    data = regrid_source_data(source, flight_number, data)

    # get the data and dimensions we're working on
    vardata = isa(var, String) ? data[var] : var
    if !isnothing(varg)
        vardatag = isa(varg, String) ? data[varg] : varg
    end
    z_dim_num = isa(z_dim, String) ? dim_num(z_dim, vardata) : z_dim # if the dim is a string, `data` must be provided so we can resolve it
    time_dim_num = isa(time_dim, String) ? dim_num(time_dim, vardata) : time_dim

    if isnothing(z_old)
        z_old = (source isa LESOutput) ? regrid_source_z_old(source, data, flight_number, interpolant_value_types.parameters[2]; thermodynamics_backend) : regrid_source_z_old(source, data; thermodynamics_backend)
    end
    if isnothing(t_old)
        if source isa LESOutput
            fix_rounding_error = true
            fix_averaging = true
            t_old = regrid_source_t_old(source, data, Val(fix_rounding_error), Val(fix_averaging)) # snap the jittery Float32 day axis to its exact 300 s grid, and extend to 0 to end
            # regrid data in time to 0 to end, so it's independent of the bc
            if fix_averaging
                lt_old = length(t_old) - 1
                # LES Data is returned as a mean over a 300s window, <150, 450, 750...>, so we extend to `deaverage` back to <0, 300, ...>, treated separately from the chosen bc
                var_data_new = Array{eltype(vardata)}(undef, ntuple(i -> i == time_dim_num ? lt_old + 1 : size(vardata, i), ndims(vardata))...)
                selectdim(var_data_new, time_dim_num, 2:lt_old) .= (selectdim(vardata, time_dim_num, 1:size(vardata, time_dim_num)-1) .+ selectdim(vardata, time_dim_num, 2:size(vardata, time_dim_num))) ./ 2
                selectdim(var_data_new, time_dim_num, 1) .= selectdim(vardata, time_dim_num, 1) - (selectdim(vardata, time_dim_num, 2) - selectdim(vardata, time_dim_num, 1)) / 2
                selectdim(var_data_new, time_dim_num, lt_old + 1) .= selectdim(vardata, time_dim_num, size(vardata, time_dim_num)) + (selectdim(vardata, time_dim_num, size(vardata, time_dim_num)) - selectdim(vardata, time_dim_num, size(vardata, time_dim_num) - 1)) / 2
                vardata = var_data_new

                if !isnothing(varg)
                    vardatag_new = Array{eltype(vardatag)}(undef, ntuple(i -> i == time_dim_num ? lt_old + 1 : size(vardatag, i), ndims(vardatag))...)
                    selectdim(vardatag_new, time_dim_num, 2:lt_old) .= (selectdim(vardatag, time_dim_num, 1:size(vardatag, time_dim_num)-1) .+ selectdim(vardatag, time_dim_num, 2:size(vardatag, time_dim_num))) ./ 2
                    selectdim(vardatag_new, time_dim_num, 1) .= selectdim(vardatag, time_dim_num, 1) - (selectdim(vardatag, time_dim_num, 2) - selectdim(vardatag, time_dim_num, 1)) / 2
                    selectdim(vardatag_new, time_dim_num, lt_old + 1) .= selectdim(vardatag, time_dim_num, size(vardatag, time_dim_num)) + (selectdim(vardatag, time_dim_num, size(vardatag, time_dim_num)) - selectdim(vardatag, time_dim_num, size(vardatag, time_dim_num) - 1)) / 2
                    vardatag = vardatag_new
                end
            end
        else
            t_old = regrid_source_t_old(source, data)
        end
    end
    initial_ind = regrid_source_initial_ind(source, data, flight_number, t_old)

    # materialize a lazy labeled variable (breaks downstream ops if left unmaterialized)
    if hasmethod(NC.dimnames, Tuple{typeof(vardata)})
        vardata = Array(vardata)
    end

    # Currently interpolates to the nearest available times, not the exact requested time; changing that would need to happen before building the vertical splines.

    if !isnothing(varg)
        assume_monotonic || error(
            "regrid_to_z_and_time: `varg` was supplied with assume_monotonic = false, which would need the " *
            "insertion index derived per column from the data; that path is not implemented. Pass " *
            "assume_monotonic = true to append the ground value at the bottom of the z axis, or insert the " *
            "ground data yourself with `combine_air_and_ground_data` and pass the combined field as `var`.",
        )
        vardata = combine_air_and_ground_data(vardata, vardatag, z_dim_num; insert_location = ground_indices) # append ground data with 0 as the z bottom; drops labeling. Assumes the two variables share size/labels, which isn't guaranteed.
    end

    # select only the initial condition timestep (keeping the dim via `[]`), else truncate from the
    # initial condition to the end along the time dimension
    if initial_condition
        vardata = selectdim(vardata, time_dim_num, [initial_ind])
    else
        vardata = selectdim(vardata, time_dim_num, initial_ind:length(t_old))
    end
    # `z_old` follows the same time selection only when it is time-varying (Atlas); LES z is static
    if regrid_source_z_time_varying(source)
        z_old =
            initial_condition ? selectdim(z_old, time_dim_num, [initial_ind]) :
            selectdim(z_old, time_dim_num, initial_ind:length(t_old))
    end

    # reverse z to run ground -> top (only for sources stored top -> ground, e.g. Atlas)
    if regrid_source_reverse_z(source)
        z_old = reverse(z_old; dims = z_dim_num) # not reverse!() bc goes SubArray --> Array
        vardata = reverse(vardata; dims = z_dim_num)
    end

    if return_before_interp
        return vardata
    end

    # evaluate onto z_new. This stage is transient — nothing is returned from it — so the coordinate
    # takes the VALUE storage spec rather than the stored-time-axis one.
    vardata = var_to_new_coord(
        vardata,
        z_old,
        z_dim_num,
        interpolant_value_types,   # transient z-interpolation: coord (irregular z) + value both use the value storage
        interpolant_value_types;
        coord_new = z_new,
        data = data,
        method = z_regrid_opts.method,
        conservative = z_regrid_opts.conservative,
        bc = z_regrid_opts.bc,
        weight = weight,
        weight_regridded = weight_regridded,
        A = A,
        Af = Af,
    )

    # caller wants just the z-regridded field (e.g. an unweighted regrid of ρ = the reusable weight denominator)
    if return_after_z_interp
        return vardata
    end

    if initial_condition # initial-condition path: nothing further to build here
        return vardata
    end

    # create new time splines. Start the times at t = 0 so the built splines are model-clock aligned.
    vardata = var_to_new_coord(
        vardata,
        t_old[initial_ind:end] .- regrid_source_t_origin(source, data, t_old, initial_ind),
        time_dim_num,
        interpolant_coord_types,   # stored time coordinate → range (O(1) per-step eval)
        interpolant_value_types,
        t_regrid_opts.drop_collinear, # prunes the returned time-spline interpolants
        ;
        coord_new = nothing,
        data = data,
        method = t_regrid_opts.method,
        conservative = t_regrid_opts.conservative,
        bc = t_regrid_opts.bc,
    )

    return vardata
end

"""
Currently this is setup to just assume saturation w/ Tg
However, this doesn't match Atlas's simulations so maybe we'll adjust this to just be an adiabatic adjustment to the lowest datapoint we do have?
It seems that Atlas's simulations have a slight kink at the lowest level, but otherwise are ≈ constant down to the surface if that's outside the forcing range and just interpolated if it's inside the forcing range...
    - in that case, we should be able to just use interpolate_1d because its default bc is nearest outside the range
"""
function calc_qg_extrapolate_pq(pg, p, q; interp_method = Interpolation.FastLinear1DInterpolation)
    # linear in p (alternative: logarithmic, i.e. linear in z not chosen, i think linear in p is right)
    qg = Interpolation.interpolate_1d(pg, p, q, interp_method; bc = Interpolation.ExtrapolateBoundaryCondition()) # default to spline1d for extrapolation as pchip may not be the most reliable outside bounds of data.
    return qg
end

"""
Retrieve the starting index for our LES simulations (12 hours before reference file) in the forcing data
"""
function initial_index(data, flight_number::Int; t_old = nothing)
    if isnothing(t_old)
        t_old = Array(data["tsec"]) # tsec from the file; verify units and that it is offset to start at 0
    end

    bdate_data = Array(data["bdate"])
    bdate_scalar = if bdate_data isa AbstractArray
        non_missing = collect(skipmissing(vec(bdate_data)))
        isempty(non_missing) && error("bdate is missing for flight $(flight_number)")
        non_missing[1]
    else
        bdate_data
    end
    bdate_str = lpad(string(Int(round(bdate_scalar))), 6, '0')
    t_base = Dates.DateTime(bdate_str, Dates.DateFormat("yymmdd")) + Dates.Year(2000) # the base Date (using bdate not nbdate cause nbdate seems to have  bug in flight 9 (extra 0 in month spot))
    t = t_base .+ Dates.Second.(t_old) # the actual dates

    initial_time = get_socrates_initial_time(flight_number, Val(false))
    initial_ind = argmin(abs.((t .- initial_time))) # find the index of the initial time
    return initial_ind
end
