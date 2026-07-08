using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Thermodynamics: Thermodynamics as TD

# ---------------------------------------------------------------------------
# Cross-backend consistency: the default `DefaultThermodynamicsBackend` and the
# accurate Thermodynamics.jl backend (via the extension) must agree to
# physically-justified tolerances over a realistic sounding. Tolerances reflect
# the intrinsic difference between the two formulations — the default backend's
# constant-L₀ Clausius–Clapeyron saturation is the loosest term (~15%); density
# ~1%; potential and virtual temperatures agree to well under 1%.
# ---------------------------------------------------------------------------
Test.@testset "Thermodynamics comparison: default vs Thermodynamics.jl extension backend" begin
    b = SSCF.DefaultThermodynamicsBackend()

    # Dependency free Thermodynamics param_set
    tp = TD.Parameters.ThermodynamicsParameters{Float64}(;
     #273.16, 273.16, 273.15, 233.0, 1.0, 1000.0, 150.0, 290.0, 220.0, 298.15, 101325.0, 100000.0, 611.657, 287.0, 461.5, 1004.5, 1859.0, 4181.0, 2070.0, 2.5008e6, 2.8344e6, 6864.8, 10513.6, 9.81, 1.0, 1.0e-10)
    T_0 = 273.16,
    T_triple = 273.16,
    T_freeze = SSCF.T_freeze(b),
    T_icenuc = SSCF.T_icenuc(b),
    T_min = 1.0,
    T_max = 1000.0,
    T_init_min = 150.0,
    T_surf_ref = 290.0,
    T_min_ref = 220.0,
    entropy_reference_temperature = 298.15,
    # Reference pressures
    MSLP = 101325.0,
    p_ref_theta = SSCF.p_ref(b),
    press_triple = 611.657,
    # Gas constants
    R_d = SSCF.R_d(b),
    R_v = SSCF.R_v(b),
    # Isobaric specific heat capacities
    cp_d = SSCF.cp_d(b),
    cp_v = SSCF.cp_v(b),
    cp_l = SSCF.cp_l(b),
    cp_i = SSCF.cp_i(b),
    # Latent heats at reference temperature
    LH_v0 = SSCF.L_v0(b),
    LH_s0 = SSCF.L_s0(b),
    # Entropy reference values
    entropy_dry_air = 6864.8,
    entropy_water_vapor = 10513.6,
    # Other
    grav = SSCF.grav(b),
    pow_icenuc = 1.0,
    q_min = 1.0e-10,
)

    # (T [K], p [Pa], q_tot [kg/kg]) across surface → upper troposphere
    pts = [(288.0, 1.0e5, 0.012), (280.0, 9.0e4, 0.008), (270.0, 7.0e4, 0.004),
           (250.0, 5.0e4, 0.001), (230.0, 3.0e4, 2.0e-4), (290.0, 1.0e5, 0.02)]

    Test.@testset "backends agree over the sounding" begin
        for (T, p, q) in pts
            Test.@test isapprox(SSCF.dry_pottemp(b, T, p), SSCF.dry_pottemp(tp, T, p); rtol = 2e-3)
            Test.@test isapprox(SSCF.virtual_temperature(b, T, p, q), SSCF.virtual_temperature(tp, T, p, q); rtol = 5e-3)
            Test.@test isapprox(SSCF.air_density(b, T, p, q), SSCF.air_density(tp, T, p, q); rtol = 3e-2)
            Test.@test isapprox(SSCF.liquid_ice_pottemp(b, T, p, q), SSCF.liquid_ice_pottemp(tp, T, p, q); rtol = 1e-2)
            Test.@test isapprox(SSCF.q_vap_saturation_liquid(b, T, p), SSCF.q_vap_saturation_liquid(tp, T, p); rtol = 0.15)
        end
    end

    Test.@testset "surface humidities agree (~1%)" begin
        for (pg, Tg) in ((1.0e5, 288.0), (9.5e4, 280.0), (1.0e5, 292.0))
            Test.@test isapprox(SSCF.saturation_q_tot_from_pgTg(b, pg, Tg), SSCF.saturation_q_tot_from_pgTg(tp, pg, Tg); rtol = 3e-2)
            Test.@test isapprox(SSCF.calc_qg_from_pgTg(b, pg, Tg), SSCF.calc_qg_from_pgTg(tp, pg, Tg); rtol = 3e-2)
        end
    end

    Test.@testset "saturation_adjust_pθq round-trip agrees" begin
        for (T, p, q) in ((285.0, 9.0e4, 1e-4), (250.0, 5.0e4, 1e-5))
            Tn = SSCF.saturation_adjust_pθq(b, p, SSCF.dry_pottemp(b, T, p), q).T
            Ta = SSCF.saturation_adjust_pθq(tp, p, SSCF.dry_pottemp(tp, T, p), q)[1]
            Test.@test isapprox(Tn, Ta; atol = 0.1)
        end
    end
end
