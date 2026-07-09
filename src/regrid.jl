#=
Regridding pipeline: build/evaluate per-column interpolants along a dimension (`interp_along_dim`,
`var_to_new_coord`), and the unified `regrid_to_z_and_time` that regrids a source onto a new vertical
grid and then builds time splines. Source-specific steps are generic verbs dispatched on
`AbstractRegridSource` (no `isa` branching in the shared body).

The interpolation API lives in the `Interpolation` submodule; call it via qualified `Interpolation.foo`.
=#

"""
    interp_along_dim

interpolation data

# var : the data to be interpolated
# interp_dim: the name or dimension number along which we will do interpolation

# interp_dim_in: the coordinate on which the data is currently aligned
# interp_dim_out: the coordinate on which we would like to interpolate the data to

the function will create splines along interp_dim_in...
- if we set interp_dim_out, we actually evaluate the function along interp_dim at locations interp_dim_out,otherwise, we just return our unevaluated spline functions

- data_func is a func to be applied to the raw data before it is processed, though perhaps it it most useful if interp_dim_out is unset and we are returning functions...
- vectorize_in means your input is an array and you need to loop over it too (as opposed to just being a fixed template vector)

Note the conservative_regirdder() returns actual output arrays -- it won't work if you are just creating a spline fcn

To Do : decide types
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
    interp_method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
    interp_kwargs::NamedTuple = (;),
    conservative_interp::Bool = false,
    conservative_interp_kwargs::Interpolation.DCIKT = Interpolation.default_conservative_interp_kwargs,
    A::Union{Nothing, AbstractArray} = nothing,
    Af::Union{Nothing, AbstractArray, LinearAlgebra.Factorization} = nothing, # precomputed factorization of A :: technically, A being Diagonal, or Triangular or something could lead to AbstractMatrix Af so we allow both AbstractMatrix and Factorization
) where {interpolant_coord_types <: Tuple, interpolant_value_types <: Tuple, drop_collinear}
    # would use kwargs but doesn't play nice w/ ODE solver for some reason... instead we get out the parameters we want explicitly and pass them all the time (splatting gives a union typle type object that apparently can't be handled?)
    f_enhancement_factor = get(interp_kwargs, :f_enhancement_factor, 1) # default to 1.0
    f_p_enhancement_factor = get(interp_kwargs, :f_p_enhancement_factor, 1) # default to 1.0
    if conservative_interp
        bc::Interpolation.AbstractBoundaryCondition = Interpolation.create_bc(get(interp_kwargs, :bc, Interpolation.ExtrapolateBoundaryCondition())) # default to extrapolate
    else
        bc = Interpolation.create_bc(get(interp_kwargs, :bc, Interpolation.ErrorBoundaryCondition())) # default to error
    end
    k = get(interp_kwargs, :k, 1) # default to 3



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
        conservative_interp ?
        ((xnew, xp, fp) -> Interpolation.conservative_regridder(
            xnew,
            xp,
            fp;
            method = interp_method,
            bc = bc,
            k = k,
            f_enhancement_factor = f_enhancement_factor,
            f_p_enhancement_factor = f_p_enhancement_factor,
            A = A,
            Af = Af,
            conservative_interp_kwargs...,
        )) : ((xnew, xp, fp) -> Interpolation.interpolate_1d(xnew, xp, fp, interp_method; bc = bc))
    # coordinate and value get INDEPENDENT storage specs (backing, eltype); Base.Fix1 fixes the spec TYPE, varies v.
    coord_conv = Base.Fix1(Interpolation.coerce_vector, interpolant_coord_types)
    value_conv = Base.Fix1(Interpolation.coerce_vector, interpolant_value_types)

    # mapslices to apply along timedim, see https://docs.julialang.org/en/v1/base/arrays/#Base.mapslices
    if !interp_dim_in_is_full_array
        if isnothing(interp_dim_out)
            # BUILD path: one interpolant per column, pruned per `drop_collinear` during construction
            # (no full-node interpolant materialized). `map` takes the array eltype from the results —
            # concrete when the backing type is length-invariant (Vector, …), abstract when it encodes
            # the node count (SVector{N}, …), for which coerce_to_shared_nodes is the opt-in unifier.
            # One path for any AbstractVector backing and either drop setting.
            xin = coord_conv(interp_dim_in)
            result = map(eachslice(vardata; dims = Tuple(d for d in 1:ndims(vardata) if d != interp_dim_num))) do d
                Interpolation.build_spline(interp_method, xin, value_conv(d); bc = bc, drop_collinear = drop_collinear_val)
            end
            return add_dim(result, interp_dim_num) # reinsert the size-1 interp dim to match the evaluate-path layout

        else
            out = mapslices(
                let interp_dim_out = interp_dim_out, interp_dim_in = interp_dim_in, interp_func = interp_func, coord_conv = coord_conv, value_conv = value_conv
                    d -> interp_func(interp_dim_out, coord_conv(interp_dim_in), value_conv(d))
                end,
                vardata,
                dims = (interp_dim_num,),
            ) # lambda fcn will evaluate, svector less meaningful

            return out
        end
    else # vectorize over input dim values as well as data (no support for vectorize over output dim yet)
        # stack on new catd dimension, then split apart inside the fcn call
        catd = ndims(vardata) + 1
        _input = cat(interp_dim_in, vardata; dims = catd) # although maybe the input is just a vector in which case this won't work..... we could just pass it in
        if isnothing(interp_dim_out)
            # BUILD path (stacked [coord data] layout): same unified map as the non-full-array case;
            # each slice is the 2-column matrix (d[:,1] = coord, d[:,2] = data), pruned per
            # `drop_collinear` during construction; `map` sets the eltype from the results.
            result = map(eachslice(_input; dims = Tuple(d for d in 1:ndims(_input) if d != interp_dim_num && d != catd))) do d
                Interpolation.build_spline(interp_method, coord_conv(d[:, 1]), value_conv(d[:, 2]); bc = bc, drop_collinear = drop_collinear_val)
            end
            return add_dim(result, interp_dim_num) # reinsert the size-1 interp dim (catd reduced out)
        else
            out = dropdims(
                mapslices(
                    let interp_dim_out = interp_dim_out, interp_func = interp_func, coord_conv = coord_conv, value_conv = value_conv
                        d -> interp_func(interp_dim_out, coord_conv(d[:, 1]), value_conv(d[:, 2]))
                    end,
                    _input,
                    dims = (interp_dim_num, catd),
                );
                dims = catd,
            ) # lambda fcn will evaluate

            return out
        end
    end
end

"""
    var_to_new_coord

if coord_new is nothing, then will return functions...
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
    interp_method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
    interp_kwargs::NamedTuple = (;),
    conservative_interp::Bool = false,
    conservative_interp_kwargs::Interpolation.DCIKT = Interpolation.default_conservative_interp_kwargs,
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

    if isnothing(weight) # if weight is a number, we assume it's a scalar and we don't need to do anything
        return interp_along_dim(
            vardata,
            interp_dim,
            coord_in,
            interpolant_coord_types,
            interpolant_value_types,
            Val(drop_collinear);
            interp_dim_out = coord_new,
            data = data,
            data_func = data_func,
            interp_dim_in_is_full_array = (size(coord_in) == size(vardata)),
            interp_method = interp_method,
            interp_kwargs = interp_kwargs,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs,
            A = A,
            Af = Af,
        )

    else

        num = interp_along_dim(
            vardata .* weight,
            interp_dim,
            coord_in,
            interpolant_coord_types,
            interpolant_value_types,
            Val(drop_collinear);
            interp_dim_out = coord_new,
            data = data,
            data_func = data_func,
            interp_dim_in_is_full_array = (size(coord_in) == size(vardata)),
            interp_method = interp_method,
            interp_kwargs = interp_kwargs,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs,
            A = A,
            Af = Af,
        )
        # Denominator = interp(weight). Reuse the caller-precomputed value if given (identical across every
        # field sharing this weight/grid/method); otherwise compute it here (bit-identical to before).
        den =
            isnothing(weight_regridded) ?
            interp_along_dim(
                weight,
                interp_dim,
                coord_in,
                interpolant_coord_types,
                interpolant_value_types,
                Val(drop_collinear);
                interp_dim_out = coord_new,
                data = data,
                data_func = data_func,
                interp_dim_in_is_full_array = (size(coord_in) == size(vardata)),
                interp_method = interp_method,
                interp_kwargs = interp_kwargs,
                conservative_interp = conservative_interp,
                conservative_interp_kwargs = conservative_interp_kwargs,
                A = A,
                Af = Af,
            ) : weight_regridded
        return num ./ den

    end
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
regrid_source_z_old(::LESOutput, data; thermodynamics_backend) =
    vec(Array(data["z"])) # static column, already present in the LES output

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
function regrid_source_t_old(::LESOutput, data, ::Val{fix_rounding_error} = Val(true)) where {fix_rounding_error} # Val for type stability 
    if fix_rounding_error
        n = length(data["time"])
        out = Vector{Int64}(undef, n) # could do like UInt32 but not worth extra readability loss
        t = vec(Array(data["time"]))
        @inbounds for i in eachindex(out)
            # out[i] = Int64(Int64(300) * (i - 1)) # fast
            out[i] = Int64(300 * round((t[i] - t[1]) * 86400 / 300)) # works on sliced data, etc
        end
    else
        t = vec(Array(data["time"])) # LES stores time in days
        out = @. (t - t[1]) * (24 * 3600) # -> seconds from start
    end
    return out
end

# index of the reference/initial timestep within `t_old`
regrid_source_initial_ind(::AtlasInput, data, flight_number, t_old) =
    initial_index(data, flight_number; t_old = t_old)
regrid_source_initial_ind(::LESOutput, data, flight_number, t_old) = 1 # LES output starts at the initial condition

# is `z_old` time-varying (needs per-timestep selection) rather than a single static column?
regrid_source_z_time_varying(::AtlasInput) = true
regrid_source_z_time_varying(::LESOutput) = false

# must the z axis be reversed to run ground -> top?
regrid_source_reverse_z(::AtlasInput) = true  # Atlas pressure levels come top -> ground
regrid_source_reverse_z(::LESOutput) = false  # LES output is already ground -> top

"""
    regrid_to_z_and_time(var, z_new, z_dim, time_dim, flight_number, interpolant_fieldtype = Vector; source = AtlasInput(), kwargs...)

Take data from a regrid `source`, interpolate it to `z_new`, then build time splines over the times
we have. To vectorize properly over `z_new` it should be the same shape as `vardata` (+ ground). The
`source` (default [`AtlasInput`](@ref); [`LESOutput`](@ref) for LES files) selects the few
source-specific steps via the `regrid_source_*` verbs; everything else is shared.
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
    interp_method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
    interp_kwargs::NamedTuple = (;),
    drop_collinear::Val = Val(false), # prune collinear nodes of the built (returned) time-spline interpolants; caller-overridable
    ground_indices = :end,
    conservative_interp::Bool = false,
    conservative_interp_kwargs::Interpolation.DCIKT = Interpolation.default_conservative_interp_kwargs,
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

    if conservative_interp && isnothing(A)
        A = Interpolation.conservative_mass_matrix(
            z_new;
            method = interp_method,
            bc = Interpolation.create_bc(get(interp_kwargs, :bc, Interpolation.ExtrapolateBoundaryCondition())),
            k = get(interp_kwargs, :k, 1),
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
            interp_method = interp_method,
            interp_kwargs = interp_kwargs,
            ground_indices = ground_indices,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs,
            A = A,
            Af = Af,
        )
    end

    data = regrid_source_data(source, flight_number, data)

    # get the data and dimensions we're working on
    vardata = isa(var, String) ? data[var] : var
    if ~isnothing(varg)
        vardatag = isa(varg, String) ? data[varg] : varg
    end
    z_dim_num = isa(z_dim, String) ? dim_num(z_dim, vardata) : z_dim # if the dim is a string, `data` must be provided so we can resolve it
    time_dim_num = isa(time_dim, String) ? dim_num(time_dim, vardata) : time_dim

    if isnothing(z_old)
        z_old = regrid_source_z_old(source, data; thermodynamics_backend)
    end
    if isnothing(t_old)
        if source isa LESOutput
            t_old = regrid_source_t_old(source, data,(_itp_coord_eltype(interpolant_coord_types) <: AbstractRange ? Val(true) : Val(true))) # For abstract ranges, go to fixed spacing (actually let's just do it no matter what)
        else
            t_old = regrid_source_t_old(source, data)
        end
    end
    initial_ind = regrid_source_initial_ind(source, data, flight_number, t_old)

    # materialize a lazy labeled variable (breaks downstream ops if left unmaterialized)
    if hasmethod(NC.dimnames, Tuple{typeof(vardata)})
        vardata = Array(vardata)
    end

    ### SHOULD WE INTERPOLATE TO THE EXACT TIME RATHER THAN CLOSEST TIMES? IDK... would need to be done before creating vertical splines...

    if ~isnothing(varg)
        # here we also are gonna need to check where things get inserted in case they are not in order...
        if !assume_monotonic # use data to figure out how and where to do insertions...
        # we need some way to get the local dimension from just a variable
        else
            vardata = combine_air_and_ground_data(vardata, vardatag, z_dim_num; insert_location = ground_indices) # append ground data with 0 as z bottom, loses labeling now though  (this puts a lot of faith im these 2 vars being the same size of having labels which we can't guarantee, no?)
        end
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

    # interpolate to new z (evaluate onto z_new; `interpolant_fieldtype` is the 4th positional arg)
    vardata = var_to_new_coord(
        vardata,
        z_old,
        z_dim_num,
        interpolant_value_types,   # transient z-interpolation: coord (irregular z) + value both use the value storage
        interpolant_value_types;
        coord_new = z_new,
        data = data,
        interp_method = interp_method,
        interp_kwargs = interp_kwargs,
        conservative_interp = conservative_interp,
        conservative_interp_kwargs = conservative_interp_kwargs,
        weight = weight,
        weight_regridded = weight_regridded,
        A = A,
        Af = Af,
    )

    # caller wants just the z-regridded field (e.g. an unweighted regrid of ρ = the reusable weight denominator)
    if return_after_z_interp
        return vardata
    end

    if initial_condition # no need to push further here since is init condition (maybe change later to return both?)
        return vardata
    end

    # create new time splines. Start the times at t = 0 so the built splines are model-clock aligned.
    vardata = var_to_new_coord(
        vardata,
        t_old[initial_ind:end] .- t_old[initial_ind],
        time_dim_num,
        interpolant_coord_types,   # stored time coordinate → range (O(1) per-step eval)
        interpolant_value_types,
        drop_collinear, # caller-set knob (default Val(false)); prunes the returned time-spline interpolants
        ;
        coord_new = nothing,
        data = data,
        interp_method = Interpolation.FastLinear1DInterpolation, # in time, we're gonna stick to linear interpolation for now... this one maybe can be pchip since it's all within the data bounds? idk... i was getting w=0 using pchip... possibly because
        interp_kwargs = interp_kwargs,
        conservative_interp = false, # no need for conservation in time? maybe?
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
    # not sure if this should be linear in p or logarithmic (linear in z), gonna do linear in p
    qg = Interpolation.interpolate_1d(pg, p, q, interp_method; bc = Interpolation.ExtrapolateBoundaryCondition()) # default to spline1d for extrapolation as pchip may not be the most reliable outside bounds of data.
    return qg
end

"""
Retrieve the starting index for our LES simulations (12 hours before reference file) in the forcing data
"""
function initial_index(data, flight_number::Int; t_old = nothing)
    if isnothing(t_old)
        t_old = Array(data["tsec"]) # check this unit was right in the files (may need to make sure it's subtracting out the first timestep so starts at 0) -- do we need to align this on a dimension?
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
