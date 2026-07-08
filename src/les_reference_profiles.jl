
"""
This is a function to get the pressure and density profiles for a given flight number/forcing to use as the reference profiles in TC.jl
    We will:
        Read p, t
        Calculate ρ
        return a named tuple with p_c, p_f, ρ_c, ρ_f

SAM seems to hold ρ constant with time, so we can probably just rely on that
"""
function les_reference_profiles(
    flight_number::Int;
    forcing_type::AbstractForcingType = ObsForcing(),
    new_zc::Union{Nothing, AbstractArray} = nothing,
    new_zf::Union{Nothing, AbstractArray} = nothing,
    # thermo_params,
)

    # initial conditions (input forcing file, opened; `.grid_data` is the vertical grid vector)
    inp = open_atlas_les_input(flight_number, forcing_type; open_files = true)

    return_values = (:p_c, :p_f, :ρ_c, :ρ_f)

    # if nothing is given, read in z from grid data to serve as new zf (plus 0 at sfc)
    # if only zf is given, calculate zc as midpoints between zf
    # if both zf and zc are given, use them as is (trust lol)
    # if only zc is given, throw an error since we don't want to do math to ensure we have repeated dzs to make that work

    # Setup new zc and zf based on input
    if isnothing(new_zc) && isnothing(new_zf)
        new_zf = [0; inp.grid_data] # by default, the grid we read in becomes zf the way TC.jl is set up
        new_zc = (new_zf[1:(end - 1)] .+ new_zf[2:end]) ./ 2
    elseif isnothing(new_zc) && isa(new_zf, AbstractArray)
        new_zc = (new_zf[1:(end - 1)] .+ new_zf[2:end]) ./ 2
    elseif isa(new_zc, AbstractArray) && isa(new_zf, AbstractArray)
        # use as-is
    else
        error("You must provide either new_zc or new_zf, or both, but only providing new_zc is not allowed")
    end

    LES_data = open_atlas_les_output(flight_number, forcing_type).data

    if isnothing(LES_data)
        error("No LES data found for flight $flight_number")
    end

    z = NC.nomissing(vec(Array(LES_data["z"])))

    p = NC.nomissing(vec(Array(LES_data["p"]))) .* 100.0 # Pressure variations in SAM are under a milibar, so we can use the 1D p variable rather than the 2D PRES variable
    ps = LES_data["Ps"][1] * 100.0 # surface pressure
    ρ = NC.nomissing(LES_data["RHO"][:, 1]) # use t=0 as our reference
    # extrapolate to get ρs since that's not given (not calculating from first principles probably is safer too w/ uncertainty in q)
    ρs = Interpolation.interpolate_1d(ps, reverse(p), reverse(ρ), Interpolation.FastLinear1DInterpolation, bc = Interpolation.ExtrapolateBoundaryCondition()) # switch to increasing for interpolation

    # Consider switching to SVector but probably not worth it for now since the high res grids can have 320 levels
    p = [ps; p]
    ρ = [ρs; ρ]
    z = [0; z]

    p_c = Interpolation.interpolate_1d(new_zc, z, p, Interpolation.FastLinear1DInterpolation; bc = Interpolation.NearestBoundaryCondition())
    p_f = Interpolation.interpolate_1d(new_zf, z, p, Interpolation.FastLinear1DInterpolation; bc = Interpolation.NearestBoundaryCondition())
    ρ_c = Interpolation.interpolate_1d(new_zc, z, ρ, Interpolation.FastLinear1DInterpolation; bc = Interpolation.NearestBoundaryCondition())
    ρ_f = Interpolation.interpolate_1d(new_zf, z, ρ, Interpolation.FastLinear1DInterpolation; bc = Interpolation.NearestBoundaryCondition())

    return NamedTuple{return_values}((p_c, p_f, ρ_c, ρ_f))

end
