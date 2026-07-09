# --- shared surface setup ------------------------------------------------------------------- #
# Load the forcing dataset and the reference-time / SST-offset metadata common to both surface
# entry points. Returns `(; data, z_dim_num, time_dim_num, initial_ind, Tg_offset)`.
function _surface_setup(flight_number::Integer, forcing_type::AbstractForcingType)
    data = open_atlas_les_input(flight_number, forcing_type).data
    z_dim_num = dim_num("lev", data["T"])
    time_dim_num = dim_num("time", data["T"])
    initial_ind = initial_index(data, flight_number)

    # deltaT is T_2m - SST (Atlas §3; table 2 has the sign backwards)
    Tg_offset = get_Tg_offset(flight_number)
    return (; data, z_dim_num, time_dim_num, initial_ind, Tg_offset)
end

wt_to_qt(wt::FT) where {FT} = wt / (one(FT) + wt) # mixing ratio to specific humidity
wt_to_qt(wt::AbstractArray{FT}) where {FT} = wt ./ (one(FT) .+ wt) # mixing ratio to specific humidity
function wt_to_qt!(wt::AbstractArray{FT}) where FT
    @. wt /= (one(FT) + wt)
    return wt
end        

# Worker: compute the reference state from an already-loaded `setup`, using `q_tot` as scratch for
# the ambient-RH branch. Both public entry points load the dataset exactly once, then call this.
function _surface_reference_state!(q_tot::AbstractArray{FT}, setup, thermodynamics_backend) where {FT}
    (; data, z_dim_num, time_dim_num, initial_ind, Tg_offset) = setup
    Tg_orig = vec(Array(data["Tg"]))[initial_ind] # SST
    Tg = Tg_orig + Tg_offset
    pg = vec(Array(data["Ps"]))[initial_ind]

    if Tg_offset < 0  # SST > Tg: SST sets q_tot at ground level and serves as a source
        q_tot_g = calc_qg_from_pgTg(thermodynamics_backend, pg, Tg)
    else # SST < Tg: stable boundary layer, ambient RH controls q_tot (not SST)
        p = vec(Array(data["lev"]))
        q_tot .= NCDatasets.nomissing(read_profile_at_time(Array(data["q"]), z_dim_num, time_dim_num, initial_ind))
        wt_to_qt!(q_tot) # mixing ratio -> specific humidity, in place
        q_tot_g = calc_qg_extrapolate_pq(pg, p, q_tot)
    end
    return (; pg = FT(pg), Tg = FT(Tg), q_tot_g = FT(q_tot_g))
end

"""
    get_SSCF_surface_reference_state!(q_tot, flight_number, forcing_type; thermodynamics_backend)

In-place surface reference state at the reference timestep, `(; pg, Tg, q_tot_g)` scalars.
`q_tot` is caller-supplied scratch for the vertical profile buffer.
"""
get_SSCF_surface_reference_state!(
    q_tot::AbstractArray,
    flight_number::Integer,
    forcing_type::AbstractForcingType;
    thermodynamics_backend = DefaultThermodynamicsBackend(),
) = _surface_reference_state!(q_tot, _surface_setup(flight_number, forcing_type), thermodynamics_backend)

"""
    get_surface_reference_state(flight_number, forcing_type, FT = Float64; thermodynamics_backend)

Allocating surface reference state at the reference timestep: `(; pg, Tg, q_tot_g)` scalars.
Loads the dataset once and forwards to [`get_SSCF_surface_reference_state!`](@ref).
"""
function get_surface_reference_state(
    flight_number::Integer,
    forcing_type::AbstractForcingType,
    ::Type{FT} = Float64;
    thermodynamics_backend = DefaultThermodynamicsBackend(),
) where {FT}
    setup = _surface_setup(flight_number, forcing_type)
    q_tot = Vector{FT}(undef, length(setup.data["lev"]))
    return _surface_reference_state!(q_tot, setup, thermodynamics_backend)
end




"""
    get_surface_forcing(flight_number, forcing_type,
                           interpolant_coord_types = Tuple{StepRangeLen, Nothing},
                           interpolant_value_types = Tuple{Vector, Float64};
                           thermodynamics_backend)

Time-dependent surface conditions from the reference timestep onward. Returns
`(; pg, Tg, Tsfc, qg, qsfc)` — each field a built time interpolant (extrapolating BC).

The `Tuple{Backing, Eltype}` storage specs (as in [`get_column_forcing`](@ref)) are threaded through
[`Interpolation.coerce_vector`](@ref) so the returned interpolants are type-stable and allocation-free
to evaluate; the default coordinate spec stores the shared time axis as a `UniformRange` (O(1) lookup).
`drop_collinear` prunes collinear nodes of the built interpolants; it must stay `false` with a
range-backed coordinate spec (pruning would break the axis's uniformity).
"""
function get_surface_forcing(
    flight_number::Integer,
    forcing_type::AbstractForcingType,
    ::Type{interpolant_coord_types} = Tuple{StepRangeLen, Nothing},
    ::Type{interpolant_value_types} = Tuple{Vector, Float64};
    thermodynamics_backend = DefaultThermodynamicsBackend(),
    drop_collinear::Val = Val(false),
) where {interpolant_coord_types <: Tuple, interpolant_value_types <: Tuple}
    (; data, z_dim_num, time_dim_num, initial_ind, Tg_offset) = _surface_setup(flight_number, forcing_type)

    pg = NCDatasets.nomissing(vec(Array(data["Ps"])))[initial_ind:end]
    Tg_orig = NCDatasets.nomissing(vec(Array(data["Tg"])))[initial_ind:end] # SST
    Tg = Tg_orig .+ Tg_offset

    if Tg_offset < 0  # SST > Tg: SST sets qg going forward and serves as a source (use full Tg)
        qg = calc_qg_from_pgTg.(thermodynamics_backend, pg, Tg)
    else # SST < Tg: stable BL, hold qg at the ambient-RH value (not the higher SST-saturation value)
        p = vec(Array(data["lev"]))
        q = NCDatasets.nomissing(
            read_profiles_over_time(Array(data["q"]), z_dim_num, time_dim_num; time_indices = initial_ind:size(Array(data["q"]), time_dim_num)),
        )
        q = collect(eachcol(q)) # list of vertical profiles, one per timestep
        q = map(mr -> mr ./ (one(eltype(mr)) .+ mr), q) # mixing ratio to specific humidity
        qg = map((pg_t, q_t) -> calc_qg_extrapolate_pq(pg_t, p, q_t), pg, q)
    end

    qg_orig = calc_qg_from_pgTg.(thermodynamics_backend, pg, Tg_orig) # q* at the actual SST

    tg = vec(Array(data["tsec"]))[initial_ind:end]
    tg = tg .- tg[1] # start at t = 0 so the interpolants are model-clock aligned

    # Coerce the shared time axis and each value series into the requested storage (backing, eltype) before
    # building — so the returned interpolants are type-stable and allocation-free to evaluate. The default
    # coordinate spec stores `tg` as a `UniformRange`, giving O(1) lookups instead of a binary search.
    coord_conv = Base.Fix1(Interpolation.coerce_vector, interpolant_coord_types)
    value_conv = Base.Fix1(Interpolation.coerce_vector, interpolant_value_types)
    tg = coord_conv(tg)
    bc = Interpolation.ExtrapolateBoundaryCondition()
    return (;
        pg = Interpolation.build_spline(Interpolation.FastLinear1DInterpolation, tg, value_conv(pg); bc, drop_collinear),
        Tg = Interpolation.build_spline(Interpolation.FastLinear1DInterpolation, tg, value_conv(Tg); bc, drop_collinear),
        Tsfc = Interpolation.build_spline(Interpolation.FastLinear1DInterpolation, tg, value_conv(Tg_orig); bc, drop_collinear), # SST, for correct sensible/latent fluxes
        qg = Interpolation.build_spline(Interpolation.FastLinear1DInterpolation, tg, value_conv(qg); bc, drop_collinear),
        qsfc = Interpolation.build_spline(Interpolation.FastLinear1DInterpolation, tg, value_conv(qg_orig); bc, drop_collinear), # q* of the actual SST
    )
end



# --- selectable forcing outputs (extensible via `Val{:symbol}` dispatch, no per-field types) ---- #
# Every field `get_column_forcing` can produce; a caller's `forcing_variables` must be a subset.
# Add an output by listing its `:symbol` here and defining `compute(::Val{:symbol}, base, tb)` (and,
# if it deviates from the defaults, `output_source` / `output_interp_kwargs` / `output_positive`).
# `base` is the shared column precompute NamedTuple (data, dims, T/p/q, pg/Tg/qg, ρ, ground_indices, FT).
"""
    supported_forcing_variables

Tuple of symbols [`get_column_forcing`](@ref) can produce. A caller's `forcing_variables`
argument must be a subset. See the user guide for field descriptions.
"""
const supported_forcing_variables =
    (:dTdt_hadv, :H_nudge, :T_nudge, :dqtdt_hadv, :qt_nudge, :subsidence, :u_nudge, :v_nudge, :ug_nudge, :vg_nudge, :dTdt_rad)

# data source: `:atlas_input` (shares the column precompute) vs `:les_output` (`:dTdt_rad` is read from
# the LES output file and regridded via a separate path, so it needs no column precompute).
output_source(::Val) = :atlas_input
output_source(::Val{:dTdt_rad}) = :les_output

# per-output conservative-regrid recipe (defaults + the few that deviate)
output_interp_kwargs(::Val) = (;)
output_interp_kwargs(::Val{:H_nudge}) = (; f_enhancement_factor = 5, f_p_enhancement_factor = 8) # keep sharp inversions (not too high -> cusps)
output_interp_kwargs(::Val{:qt_nudge}) = (; f_enhancement_factor = 6, f_p_enhancement_factor = 8) # keep sharp inversions
output_interp_kwargs(::Val{:subsidence}) = (; f_enhancement_factor = 1, f_p_enhancement_factor = 1) # gentle changes; accuracy loss ok

output_positive(::Val) = false
output_positive(::Val{:H_nudge}) = true # θ_liq_ice is positive-definite
output_positive(::Val{:qt_nudge}) = true # total specific humidity is positive-definite

# pre-regrid, full-grid (air + ground) field for an output, from the shared `base` column precompute.
_full(base, air, ground) = combine_air_and_ground_data(air, ground, base.z_dim_num; insert_location = base.ground_indices)
compute(::Val{:dTdt_hadv}, base, tb) = _full(base, base.data["divT"], zero(base.FT)) # horizontal advective T tendency (ERA5)
compute(::Val{:dqtdt_hadv}, base, tb) = _full(base, base.data["divq"], zero(base.FT)) # horizontal advective qt tendency (ERA5)
compute(::Val{:u_nudge}, base, tb) = _full(base, base.data["u"], zero(base.FT)) # ERA5 horizontal wind (ground value 0)
compute(::Val{:v_nudge}, base, tb) = _full(base, base.data["v"], zero(base.FT))
compute(::Val{:ug_nudge}, base, tb) = _full(base, base.data["ug"], zero(base.FT)) # ERA5 geostrophic wind
compute(::Val{:vg_nudge}, base, tb) = _full(base, base.data["vg"], zero(base.FT))
compute(::Val{:qt_nudge}, base, tb) = _full(base, base.q, base.qg) # total specific humidity
compute(::Val{:T_nudge}, base, tb) = _full(base, base.T, base.Tg) # absolute temperature
compute(::Val{:H_nudge}, base, tb) =
    liquid_ice_pottemp.(tb, _full(base, base.T, base.Tg), _full(base, base.p, base.pg), _full(base, base.q, base.qg))
compute(::Val{:subsidence}, base, tb) =
    _column_subsidence(base.data, base.ρ, _full(base, base.p, base.pg), base.z_dim_num, base.ground_indices, tb, base.FT)





# Working float type from a value storage spec `Tuple{Backing, Eltype}`: the eltype, or `Float64` when the
# eltype knob is `Nothing` (as-read; NetCDF reads Float64). Drives thermo/A-matrices/materialized arrays.
_itp_value_eltype(::Type{<:Tuple{Any, Nothing}}) = Float64 
_itp_value_eltype(::Type{Tuple{CT,VT}}) where {CT,VT} = VT

_itp_coord_eltype(::Type{Tuple{CT,VT}}) where {CT,VT} = CT


"""
    get_column_forcing(
        flight_number,
        forcing_type,
        forcing_variables = supported_forcing_variables,
        interpolant_coord_types = Tuple{StepRangeLen, Nothing},
        interpolant_value_types = Tuple{Vector, Float64},
        FT = _itp_value_eltype(interpolant_value_types);
        new_z = nothing,
        initial_condition = false,
        thermodynamics_backend = DefaultThermodynamicsBackend(),
        use_LES_output_for_z = false,
        return_old_z = false,
        fail_on_missing_data = true,
        conservative_interp = false,
        conservative_interp_kwargs = Interpolation.default_conservative_interp_kwargs,
        drop_collinear = Val(false),
        A_cache = Dict{DataType, Matrix{FT}}(),
        Af_cache = Dict{DataType, LinearAlgebra.Factorization{FT}}(),
    )

Build column forcing for a SOCRATES flight and [`AbstractForcingType`](@ref) (`ObsForcing()` / `ERA5Forcing()`).

Returns a `NamedTuple` keyed by `forcing_variables` (default [`supported_forcing_variables`](@ref)).
Each value is a `Vector` of time interpolants (one per vertical level on `new_z`), unless
`initial_condition = true` (time-0 profiles) or `return_old_z = true` (source altitude field).

Storage type parameters `interpolant_coord_types` and `interpolant_value_types` are
`Tuple{Backing, Eltype}` specs passed to [`Interpolation.coerce_vector`](@ref); the value
spec eltype sets the working float type `FT`.

Winds, subsidence, and geostrophic fields are always ERA5-sourced; nudging targets and
advective tendencies follow `forcing_type`. Equilibrium-derived quantities use
`thermodynamics_backend` (default [`DefaultThermodynamicsBackend`](@ref); pass
`ThermodynamicsParameters` with `Thermodynamics.jl` loaded for accurate physics).
"""
function get_column_forcing(
    flight_number::FNT,
    forcing_type::AbstractForcingType,
    forcing_variables::NTuple{NFV, Symbol} = supported_forcing_variables,
    ::Type{interpolant_coord_types} = Tuple{StepRangeLen, Nothing}, # stored time-COORD spec (backing, eltype): range → O(1) eval; eltype as-read
    ::Type{interpolant_value_types} = Tuple{Vector, Float64},       # stored VALUE spec (backing, eltype); its eltype is the working type FT
    ::Type{FT} = _itp_value_eltype(interpolant_value_types),
    ;
    new_z::Union{Nothing, AbstractArray, NamedTuple} = nothing,
    initial_condition::Bool = false,
    thermodynamics_backend = DefaultThermodynamicsBackend(),
    use_LES_output_for_z::Bool = false,
    return_old_z::Bool = false,
    fail_on_missing_data::Bool = true,
    conservative_interp::Bool = false, # off by default
    conservative_interp_kwargs::Interpolation.DCIKT = Interpolation.default_conservative_interp_kwargs,
    # Prune collinear nodes of the built (returned) time-spline interpolants. Off by default: pruning is
    # an exact (value-preserving) size optimization but yields per-length node sets, so pass `Val(true)`
    # only with a length-invariant backing (e.g. `Vector`) where the result stays concrete.
    drop_collinear::Val = Val(false),
    # Reusable (across calls) conservative-interp caches, concretely typed: keyed by interpolation
    # *method type* (concrete `DataType`), holding the FT mass matrix and its LU factorization. `A`'s
    # eltype follows the target grid; stored as `Matrix{FT}` (a no-op when the grid is already `FT`).
    A_cache::Dict{DataType, Matrix{FT}} = Dict{DataType, Matrix{FT}}(),
    # `factorize` (kept for its adaptivity — e.g. Bunch-Kaufman on a symmetric mass matrix) can return
    # different factorization kinds, so the value is typed `Factorization{FT}`: concrete in precision,
    # abstract only in kind (one `\` dispatch per field at solve time — negligible).
    Af_cache::Dict{DataType, LinearAlgebra.Factorization{FT}} = Dict{DataType, LinearAlgebra.Factorization{FT}}(),
) where {FNT <: Integer, interpolant_coord_types <: Tuple, interpolant_value_types <: Tuple, FT, NFV}
    # working float type = the value storage's eltype (thermo, A matrices, materialized arrays). `Nothing`
    # (as-read) resolves to Float64 (the NetCDF read type).
    # `enforce_positivity` is a property of the conservative interpolation (not the forcing_type). Only
    # `qt_nudge` and `H_nudge` are positive-definite, so we default it off here and re-enable it for
    # just those two below.
    enforce_positivity = conservative_interp_kwargs[:enforce_positivity]
    conservative_interp_kwargs = set_property(conservative_interp_kwargs, :enforce_positivity, false)

    # `forcing_variables` is the single source of truth: it selects which per-field grids are built,
    # which fields are actually computed, and the returned NamedTuple. Reject anything we can't produce.
    for v in forcing_variables
        v in supported_forcing_variables ||
            error("unsupported forcing variable :$v — supported: $(supported_forcing_variables)")
    end

    data = _open_atlas_les_input(flight_number, forcing_type, Val(true), Val(new_z isa Nothing)) # load only the forcing source we need

    # broadcast the requested grid to a per-field NamedTuple keyed by `forcing_variables`; every branch
    # yields a NamedTuple{forcing_variables} so downstream indexing (`new_z[:field]`) stays consistent.
    new_z = if new_z isa NamedTuple
        new_z
    elseif new_z isa Nothing # better than isnothing() for JET.jl
        NamedTuple{forcing_variables}(ntuple(_ -> data[:grid_data], Val(length(forcing_variables))))
    else # AbstractArray
        NamedTuple{forcing_variables}(ntuple(_ -> new_z, Val(length(forcing_variables))))
    end

    data = data.data # keep only the opened dataset for this forcing source

    # specify our action dimensions
    z_dim_num = dim_num("lev", data["T"])
    time_dim_num = dim_num("time", data["T"])

    # === conservative mass-response matrix cache ================================================== #
    # A (and its factorization) are expensive integrals; compute once per interpolation method and
    # reuse across fields (and across calls, via the passed-in caches). `_AAf` is a no-op returning
    # `(nothing, nothing)` unless `conservative_interp` is set; otherwise it self-seeds from the
    # requesting field's grid. Keyed by method *type* (concrete `DataType`), and it returns the values
    # read back from the caches so they carry the concrete cache eltype (`Matrix{FT}` / its `LU`).
    function _AAf(A_cache, Af_cache, name, method)
        conservative_interp || return (nothing, nothing)
        key = typeof(method)
        if !haskey(A_cache, key)
            A_cache[key] = Interpolation.conservative_mass_matrix(
                new_z[name];
                method = method,
                k = 1,
                bc = Interpolation.ExtrapolateBoundaryCondition(),
            )
            Af_cache[key] = LinearAlgebra.factorize(A_cache[key])
        end
        return (A_cache[key], Af_cache[key])
    end

    computed = (;) # accumulate only the requested fields (concrete-typed), reordered at return

    # === Atlas-input fields (share the column precompute) ========================================= #
    # Everything except the LES-sourced `:dTdt_rad`. `return_old_z` also needs this precompute (z_old).
    needs_column = return_old_z || any(v -> output_source(Val(v)) === :atlas_input, forcing_variables)
    if needs_column
        # p,T,q — materialize once into Array{FT} so downstream broadcasts fuse over plain arrays
        p = align_along_dimension(_materialize(data["lev"]), z_dim_num)
        pg = _materialize(data["Ps"])
        T = _materialize(data["T"]) # NOTE: Atlas stores the liquid-ice temperature T_L here, not the actual T
        Tg = _materialize(data["Tg"])
        q = wt_to_qt!(_materialize(data["q"])) # mixing ratio -> specific humidity, in place on the owned copy

        # Recover the ACTUAL state from the Atlas liquid-ice temperature `T_L`: total water is conserved,
        # so saturation-adjust `(p, θ_liq_ice ≈ dry_pottemp(T_L), q_tot)` to the equilibrium
        # `(T, q_liq, q_ice)`. This reproduces the original pipeline's phase-equilibrium state; without it
        # the altitude integration uses `T_L` for `T` (biasing `z` at saturated levels) and density
        # ignores the condensate loading (`q_vap = q_tot - q_liq - q_ice`). `q_tot` is unchanged.
        _state = saturation_adjust_pθq.(thermodynamics_backend, p, dry_pottemp.(thermodynamics_backend, T, p), q)
        T = map(s -> s.T, _state)      # actual temperature
        q_liq = map(s -> s.q_liq, _state)  # equilibrium liquid (condensate-aware density)
        q_ice = map(s -> s.q_ice, _state)  # equilibrium ice

        # add the ΔT from the summary table (Atlas §3: really T_2m - SST; table 2 has the sign backwards)
        Tg_orig = Tg # SST; Tg is reassigned to a fresh array below (never mutated), so an alias is safe
        Tg_offset = get_Tg_offset(flight_number)

        @. Tg += Tg_offset

        if Tg_offset < 0  # SST > Tg, assume Tg limits moisture below last known point and just extrapolate
            qg = calc_qg_extrapolate_pq.(
                pg,
                Ref(vec(Array(p))),
                align_along_dimension(vec.(collect(eachslice(Array(q); dims = time_dim_num))), z_dim_num),
            ) # pg is a scalar-per-column, p is a fixed profile, q slices in z are aligned along time to match pg's shape
        else # Tg > SST, assume SST sets moisture below last known point
            qg = calc_qg_from_pgTg.(thermodynamics_backend, pg, Tg_orig) # SST (Tg_orig) sets qg below the last known point
        end

        p_grid = add_dim(align_along_dimension(p, z_dim_num), time_dim_num)

        # `p_grid` (lon,lat,lev,time) — NOT the 3-D `align_along_dimension(p, ...)` — so `pg` (lon,lat,time)
        # is one ndim short and gets a singleton *level* dim inserted; otherwise its time axis is
        # mistaken for extra levels and concatenated with them.
        ground_indices = ground_insertion_indices(
            p_grid,
            pg,
            z_dim_num;
            data = data,
        )

        # old_z: precompute once (shared by every column field; regrid_to_z_and_time could self-derive it
        # but it's redundant to recompute per field)
        if !use_LES_output_for_z
            z_old = lev_to_z(
                p_grid,
                T,
                q,
                pg,
                Tg,
                qg;
                data = data,
                thermodynamics_backend = thermodynamics_backend,
            )
        else
            z_old = lev_to_z_from_LES_output(
                p_grid,
                T,
                q,
                pg,
                Tg,
                qg;
                data = data,
                thermodynamics_backend = thermodynamics_backend,
                flight_number = flight_number,
                forcing_type = forcing_type,
                ground_indices = ground_indices,
            )
        end

        if return_old_z
            return z_old
        end

        # density on the full (air+ground) grid — the interpolation weight for every column field, and
        # the subsidence denominator. Uses the condensate-aware 6-arg form so `ρ` reflects the vapor
        # fraction `q_tot - q_liq - q_ice` (matters at saturated levels; the surface partition is resolved
        # from its own equilibrium condensate).
        ρ = air_density.(thermodynamics_backend, T, p, q, q_liq, q_ice)
        _cg = equilibrium_condensate.(thermodynamics_backend, Tg, pg, qg)
        ρg = air_density.(thermodynamics_backend, Tg, pg, qg, getproperty.(_cg, :q_liq), getproperty.(_cg, :q_ice))
        ρ = combine_air_and_ground_data(ρ, ρg, z_dim_num; insert_location = ground_indices)

        # The mass-weight normalization denominator `interp(ρ)` is identical for every field regridded onto
        # the same grid with the same method, so compute it ONCE here (an unweighted regrid of ρ, returned
        # right after the z-interpolation) and reuse it across all fields via `weight_regridded`, instead of
        # re-regridding ρ per field. It is field-independent only in the non-conservative path (conservative
        # regridding folds in per-field `f_enhancement`), and only when every atlas field shares one grid;
        # otherwise `nothing` → each field computes its own denominator (bit-identical to no hoisting). This
        # is a typed local (not an `Any` cache): ρ is data that changes every call, so there is nothing to
        # reuse across calls (unlike the grid-only `A_cache`/`Af_cache`).
        _atlas_vars = filter(v -> output_source(Val(v)) === :atlas_input, forcing_variables)
        weight_regridded_shared =
            (!conservative_interp && !isempty(_atlas_vars) &&
             all(v -> new_z[v] === new_z[first(_atlas_vars)], _atlas_vars)) ?
            regrid_to_z_and_time(
                ρ, new_z[first(_atlas_vars)], z_dim_num, time_dim_num, flight_number,
                interpolant_coord_types, interpolant_value_types;
                z_old = z_old, data = data, thermodynamics_backend, initial_condition,
                weight = nothing, conservative_interp = false,
                conservative_interp_kwargs = conservative_interp_kwargs,
                return_after_z_interp = true,
            ) : nothing

        # regrid one already-assembled column field onto its requested grid (+ time splines if not
        # initial_condition). By default we conservatively interpolate and weight by density — slightly
        # diffusive even for intensive variables, but it rounds corners and conserves mass.
        function _regrid(
            field,
            name;
            interp_method = Interpolation.FastLinear1DInterpolation,
            interp_kwargs = (;),
            cons_kwargs = conservative_interp_kwargs,
            A_cache = A_cache,
            Af_cache = Af_cache,
        )
            A, Af = _AAf(A_cache, Af_cache, name, interp_method)
            return regrid_to_z_and_time(
                field,
                new_z[name],
                z_dim_num,
                time_dim_num,
                flight_number,
                interpolant_coord_types,
                interpolant_value_types;
                z_old = z_old,
                data = data,
                thermodynamics_backend,
                initial_condition,
                interp_method = interp_method,
                interp_kwargs = interp_kwargs,
                drop_collinear = drop_collinear,
                conservative_interp = conservative_interp,
                conservative_interp_kwargs = cons_kwargs,
                weight = ρ,
                weight_regridded = weight_regridded_shared,
                A = A,
                Af = Af,
            )[:]
        end

        # `base.p` is the 4-D `p_grid` (lon,lat,lev,time) so full-column fields (H_nudge, subsidence)
        # combine it with the per-time 4-D `ground_indices`; the 3-D `p` is only for the `ρ` broadcast above.
        base = (; data, z_dim_num, time_dim_num, T, p = p_grid, q, pg, Tg, qg, ρ, ground_indices, FT)

        # compute + regrid each requested :atlas_input output via Val-dispatch on its symbol; the
        # `compute` methods above assemble each pre-regrid field from `base`, and
        # output_interp_kwargs / output_positive give its per-field regrid recipe.
        for v in forcing_variables
            output_source(Val(v)) === :atlas_input || continue
            field = compute(Val(v), base, thermodynamics_backend)
            cons =
                output_positive(Val(v)) ?
                set_property(conservative_interp_kwargs, :enforce_positivity, enforce_positivity) :
                conservative_interp_kwargs
            computed = merge(
                computed,
                NamedTuple{(v,)}((_regrid(field, v; interp_kwargs = output_interp_kwargs(Val(v)), cons_kwargs = cons),)),
            )
        end
    end

    # === LES-sourced field: dTdt_rad ============================================================== #
    if :dTdt_rad in forcing_variables
        LES_data = open_atlas_les_output(flight_number, forcing_type).data
        # `_materialize` bulk-loads the lazy CFVariable ONCE (the package's standard NC boundary helper);
        # operating on the lazy variable directly indexes it element-by-element (per-element allocation,
        # the same trap as `combine`).
        dTdt_rad = _materialize(LES_data["RADQR"]) ./ (FT(24) * FT(3600)) # K/day -> K/s (drops nc dim metadata, yields Array)
        z_dim_num_LES = dim_num("z", LES_data["RADQR"]) # assume same for all LES 2D vars
        time_dim_num_LES = dim_num("time", LES_data["RADQR"])
        ρ_LES = _materialize(LES_data["RHO"])

        A, Af = _AAf(A_cache, Af_cache, :dTdt_rad, Interpolation.FastLinear1DInterpolation)
        dTdt_rad = regrid_to_z_and_time(
            dTdt_rad,
            new_z[:dTdt_rad],
            z_dim_num_LES,
            time_dim_num_LES,
            flight_number,
            interpolant_coord_types,
            interpolant_value_types;
            source = LESOutput(forcing_type),
            interp_kwargs = (; bc = Interpolation.ExtrapolateBoundaryCondition()), # RF09 LES ran on a different grid -> extrapolate
            z_old = nothing,
            data = LES_data,
            initial_condition,
            drop_collinear = drop_collinear,
            conservative_interp = conservative_interp,
            conservative_interp_kwargs = conservative_interp_kwargs, # can be negative
            weight = ρ_LES, # mass-weighted conservative regrid
            A = A,
            Af = Af,
        )[:]
        computed = merge(computed, (; dTdt_rad = dTdt_rad))
    end

    # reorder to `forcing_variables` (preserving each field's concrete type) and return
    return NamedTuple{forcing_variables}(computed)

end

# Subsidence (full-grid, prior to vertical regridding) from omega and the surface pressure tendency,
# using the Atlas scaled-sigmoid weighting `f_p` that passes through (ps, 1) and (25000 Pa, 0). `ρ` is
# the full-grid air density and `p_full` the full-grid pressure; both are supplied by the caller.
function _column_subsidence(
    data,
    ρ,
    p_full,
    z_dim_num,
    ground_indices,
    thermodynamics_backend,
    ::Type{FT},
) where {FT}
    # Materialize the raw NCDataset variables at the boundary: `dpdt_g` is broadcast directly below (it
    # does not pass through `combine_air_and_ground_data`'s materialize guard), so leaving it lazy pushes a
    # `CFVariable`/`ReshapedDiskArray` into the subsidence broadcast — whose eltype infers to `Any` on some
    # Julia versions, breaking downstream. Loading once here also avoids lazy per-element disk reads.
    ω = _materialize(data["omega"])
    dpdt_g = _materialize(data["Ptend"])
    dpdt_g = add_dim(dpdt_g, z_dim_num) # should be lon lat lev time (hopefully order was already correct)
    ω = combine_air_and_ground_data(ω, dpdt_g, z_dim_num; insert_location = ground_indices)

    # Despite matching documentation, this was wrong -- Atlas confirmed what version of the sigmoid she used (see comments below)
    # c    = 100
    # a    = -2 + exp(250/c)
    # f_p  = @. 2(a+1) / (a+exp(p_full/c)) - 1

    # Weighting through (ps, 1) and (25000 Pa, 0). We use the cosine form Atlas gave (by email); the
    # logistic form below is the original derivation, kept (commented) for reference.
    p0 = FT(250.0) * FT(100.0)
    p1 = maximum(p_full, dims = z_dim_num) # ground pressure per column
    f_p = @. cos(FT(π) / FT(2) * (p1 - p_full) / (p1 - p0)) * (p_full >= p0) # atlas email (used)

    # --- logistic sigmoid alternative (unused; kept for reference) --- #
    # L = FT(2.2) # maximum value (shape parameter)
    # a = -L / FT(2)
    # (p1, y1) = (p1, FT(1)) # value 1 at the ground pressure
    # (p2, y2) = (p0, FT(0)) # value 0 at p0 = 25000 Pa
    # k = @. log((L / FT(2) + FT(1)) / (L / FT(2) - FT(1))) / (p1 - p0) # logistic rate (array over columns)
    # f_p_alt = @. (a + L / (FT(1) + exp(-k * (p_full - p0)))) * (p_full >= p0) # my original

    # `f_p`, `ω`, `dpdt_g`, and `ρ` are all built on the same full `(lon, lat, lev+1, time)` grid from
    # the same forcing file, so the subsidence just broadcasts — no dim re-alignment needed. (If a
    # future source broke that assumption, this broadcast errors loudly rather than silently permuting
    # dims to a same-size-but-wrong-semantics layout.)
    g = grav(thermodynamics_backend, FT)
    subsidence = f_p
    @. subsidence = -(ω - (dpdt_g * f_p)) / (ρ * g)
    return subsidence
end

"""
    default_new_z(flight_number)

Return the Atlas default vertical grid vector for `flight_number` (from `open_atlas_les_grid`).
"""
function default_new_z(flight_number::Int;)
    data = open_atlas_les_grid(flight_number)
    new_z = data[:grid_data]
    return new_z
end


