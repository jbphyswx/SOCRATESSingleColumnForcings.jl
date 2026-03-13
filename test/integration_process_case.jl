@warn "If you're on HPC, you should set `export OPENBLAS_NUM_THREADS=1` in your environment, or use `LinearAlgebra.BLAS.set_num_threads(1)` in your code to avoid performance issues."

@testset "SOCRATESSingleColumnForcings integration" begin
    FT = Float64

    @info "SOCRATESSingleColumnForcings version: $(pkgversion(SSCF))"
    @info "Thermodynamics version: $(pkgversion(TD))"

    toml_dict = CP.create_toml_dict(FT)
    thermo_params = TDP.ThermodynamicsParameters(toml_dict)

    data = SSCF.open_atlas_les_input(9)
    new_z = data[:grid_data]
    data = data[:obs_data]

    new_zs = (nothing, FT.(1:100:4000))

    initial_conditions = (false, true)
    surfaces = (nothing, "reference_state", "surface_conditions")
    use_LES_output_for_zs = (false, true)
    return_old_zs = (false, true)
    conservative_interps = (false, true)
    enforce_positivitys = (false, true)
    use_svectors = (true, true)

    SSCF.conservative_spline_values(
        new_zs[2],
        rand(FT, length(new_zs), 10);
        bc = SSCF.ExtrapolateBoundaryCondition(),
        enforce_positivity = true,
    )

    setups = Iterators.product(
        SSCF.flight_numbers,
        SSCF.forcing_types,
        new_zs,
        initial_conditions,
        surfaces,
        use_LES_output_for_zs,
        return_old_zs,
        conservative_interps,
        enforce_positivitys,
        use_svectors,
    )

    n_setups = length(setups)
    for (i, setup) in enumerate(setups)
        flight_number,
        forcing_type,
        new_z,
        initial_condition,
        surface,
        use_LES_output_for_z,
        return_old_z,
        conservative_interp,
        enforce_positivity,
        use_svectors = setup

        @info "Testing combination $i/$n_setups: flight_number=$flight_number, forcing_type=$forcing_type, initial_condition=$initial_condition, surface=$(surface), use_LES_output_for_z=$use_LES_output_for_z, return_old_z=$return_old_z, conservative_interp=$conservative_interp, enforce_positivity=$enforce_positivity, use_svectors=$use_svectors"

        data = SSCF.open_atlas_les_input(flight_number, forcing_type; open_files = false)
        conservative_interp_kwargs = SSCF.get_conservative_interp_kwargs(;)

        if !conservative_interp && enforce_positivity
            continue
        end
        if conservative_interp
            conservative_interp_kwargs =
                SSCF.set_property(conservative_interp_kwargs, :enforce_positivity, enforce_positivity)
        end

        if (!isnothing(data[forcing_type]) && isfile(data[forcing_type])) &&
           (!isnothing(data[:grid_data]) && isfile(data[:grid_data]))
            SSCF.process_case(
                flight_number;
                forcing_type = forcing_type,
                thermo_params = thermo_params,
                new_z = new_z,
                initial_condition = initial_condition,
                surface = surface,
                use_LES_output_for_z = use_LES_output_for_z,
                return_old_z = return_old_z,
                conservative_interp = conservative_interp,
                conservative_interp_kwargs = conservative_interp_kwargs,
                fail_on_missing_data = false,
                use_svectors = use_svectors,
            )
            @test true
        else
            @warn "Files for flight $flight_number do not exist"
        end
    end

    ref = SSCF.get_LES_reference_profiles(9; forcing_type = :obs_data, new_zc = nothing, new_zf = nothing)
    @test !isnothing(ref)

    data = SSCF.process_case(9; thermo_params = thermo_params, use_LES_output_for_z = false)
    @test !isnothing(data)

    data = SSCF.process_case(9; thermo_params = thermo_params, initial_condition = true, use_LES_output_for_z = false)
    @test !isnothing(data)

    data = SSCF.process_case(
        9;
        thermo_params = thermo_params,
        surface = "reference_state",
        use_LES_output_for_z = false,
    )
    @test !isnothing(data)

    data = SSCF.process_case(
        9;
        thermo_params = thermo_params,
        surface = "surface_conditions",
        use_LES_output_for_z = false,
    )
    @test !isnothing(data)
end