
"""
    les_reference_profiles(flight_number, FT = Float64; forcing_type = ObsForcing(), new_zc = nothing, new_zf = nothing,
                           interp_method = Interpolation.FastLinear1DInterpolation, interp_kwargs = (; bc = …),
                           conservative_interp = false, conservative_interp_kwargs = …)

Pressure and density reference profiles on the cell (`new_zc`) and face (`new_zf`) vertical grids,
returned as `FT` vectors `(; p_c, p_f, ρ_c, ρ_f)`. With neither grid given, the input file's vertical
grid becomes `new_zf` (0 prepended at the surface) and `new_zc` its midpoints; with only `new_zf`,
`new_zc` is its midpoints.

The vertical regrid uses the same interpolation machinery as the forcing profiles ([`get_column_forcing`](@ref)):
`interp_method` / `interp_kwargs` (including `bc`) and optional `conservative_interp` /
`conservative_interp_kwargs`, threaded through `var_to_new_coord`. Allocates the four output
buffers and forwards to [`les_reference_profiles!`](@ref) (call that directly to reuse your own buffers).
"""
function les_reference_profiles(
    flight_number::Int,
    ::Type{FT} = Float64;
    forcing_type::AbstractForcingType = ObsForcing(),
    new_zc::Union{Nothing, AbstractArray} = nothing,
    new_zf::Union{Nothing, AbstractArray} = nothing,
    interp_method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
    interp_kwargs::NamedTuple = (; bc = Interpolation.NearestBoundaryCondition()),
    conservative_interp::Bool = false,
    conservative_interp_kwargs::Interpolation.DCIKT = Interpolation.default_conservative_interp_kwargs,
) where {FT}
    # Resolve the target cell/face grids (only the default path reads the input file's grid vector):
    #  - neither given: the read-in grid becomes the face grid (0 at the surface), the cell grid its midpoints
    #  - only zf given: zc is its midpoints
    #  - both given: used as-is  (only zc given is disallowed — would require guessing the face spacing)
    if isnothing(new_zc) && isnothing(new_zf)
        new_zf = [0; open_atlas_les_input(flight_number, forcing_type; open_files = true).grid_data]
        new_zc = (new_zf[1:(end - 1)] .+ new_zf[2:end]) ./ 2
    elseif isnothing(new_zc) && isa(new_zf, AbstractArray)
        new_zc = (new_zf[1:(end - 1)] .+ new_zf[2:end]) ./ 2
    elseif isa(new_zc, AbstractArray) && isa(new_zf, AbstractArray)
        # use as given
    else
        error("You must provide either new_zc or new_zf, or both, but only providing new_zc is not allowed")
    end

    p_c = Vector{FT}(undef, length(new_zc)); ρ_c = similar(p_c)
    p_f = Vector{FT}(undef, length(new_zf)); ρ_f = similar(p_f)
    return les_reference_profiles!(
        p_c, p_f, ρ_c, ρ_f, flight_number;
        forcing_type, new_zc, new_zf, interp_method, interp_kwargs, conservative_interp, conservative_interp_kwargs,
    )
end

"""
    les_reference_profiles!(p_c, p_f, ρ_c, ρ_f, flight_number; forcing_type = ObsForcing(), new_zc, new_zf,
                            interp_method = …, interp_kwargs = …, conservative_interp = false, conservative_interp_kwargs = …)

In-place [`les_reference_profiles`](@ref): regrid the pressure/density reference profiles onto the given
`new_zc` (cell) and `new_zf` (face) grids into the caller-provided buffers, returning `(; p_c, p_f, ρ_c, ρ_f)`.
The buffers' element type sets the output type. The regrid goes through `var_to_new_coord`, so it
honors `interp_method` / `interp_kwargs` / `conservative_interp` exactly like the forcing profiles.
"""
function les_reference_profiles!(
    p_c::AbstractVector, p_f::AbstractVector, ρ_c::AbstractVector, ρ_f::AbstractVector,
    flight_number::Int;
    forcing_type::AbstractForcingType = ObsForcing(),
    new_zc::AbstractArray,
    new_zf::AbstractArray,
    interp_method::Interpolation.AbstractInterpolationMethod = Interpolation.FastLinear1DInterpolation,
    interp_kwargs::NamedTuple = (; bc = Interpolation.NearestBoundaryCondition()),
    conservative_interp::Bool = false,
    conservative_interp_kwargs::Interpolation.DCIKT = Interpolation.default_conservative_interp_kwargs,
)
    LES_data = open_atlas_les_output(flight_number, forcing_type).data
    isnothing(LES_data) && error("No LES data found for flight $flight_number")

    z = NC.nomissing(vec(Array(LES_data["z"])))
    FT = eltype(z)
    p = NC.nomissing(vec(Array(LES_data["p"]))) .* FT(100) # SAM's sub-millibar pressure variations let us use the 1D p over the 2D PRES
    ps = LES_data["Ps"][1] * FT(100) # surface pressure
    ρ = NC.nomissing(LES_data["RHO"][:, 1]) # t = 0 reference
    # ρs is not stored; extrapolate to the surface (increasing-pressure order for the interpolant)
    ρs = Interpolation.interpolate_1d(ps, reverse(p), reverse(ρ), Interpolation.FastLinear1DInterpolation; bc = Interpolation.ExtrapolateBoundaryCondition())

    z = [0; z]; p = [ps; p]; ρ = [ρs; ρ]

    # Regrid over `z` with the SAME machinery the forcing profiles use (`var_to_new_coord` → configurable
    # `interp_method`/`bc`/conservative). Writing into the caller's buffers sets the output element type.
    rg(field, new_z) = var_to_new_coord(
        field, z, 1;
        coord_new = new_z, interp_method, interp_kwargs, conservative_interp, conservative_interp_kwargs,
    )
    p_c .= rg(p, new_zc); p_f .= rg(p, new_zf)
    ρ_c .= rg(ρ, new_zc); ρ_f .= rg(ρ, new_zf)
    return (; p_c, p_f, ρ_c, ρ_f)
end
