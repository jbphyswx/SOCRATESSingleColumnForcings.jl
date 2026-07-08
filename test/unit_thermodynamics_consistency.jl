using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Thermodynamics: Thermodynamics as TD
using ClimaParams: ClimaParams as CP

# ---------------------------------------------------------------------------
# Cross-backend consistency: the naive `DefaultThermodynamicsBackend` and the
# accurate Thermodynamics.jl backend (via the extension) must agree to
# physically-justified tolerances over a realistic sounding. Tolerances reflect
# the intrinsic difference between the two formulations — the naive backend's
# constant-L₀ Clausius–Clapeyron saturation is the loosest term (~15%); density
# ~1%; potential and virtual temperatures agree to well under 1%.
# ---------------------------------------------------------------------------
Test.@testset "Thermodynamics seam: naive vs accurate backend" begin
    b = SSCF.DefaultThermodynamicsBackend()
    tp = TD.Parameters.ThermodynamicsParameters(CP.create_toml_dict(Float64))

    Test.@testset "physical constants agree (<0.1%)" begin
        Test.@test isapprox(SSCF.R_d(b), SSCF.R_d(tp); rtol = 1e-3)
        Test.@test isapprox(SSCF.R_v(b), SSCF.R_v(tp); rtol = 1e-3)
        Test.@test isapprox(SSCF.grav(b), SSCF.grav(tp); rtol = 1e-3)
        # ε must be the meteorological ratio M_v/M_d ≈ 0.622 in both backends (Thermodynamics' own
        # `molmass_ratio` is the reciprocal M_d/M_v ≈ 1.608, which must NOT leak through the seam)
        Test.@test isapprox(SSCF.molmass_ratio(b), SSCF.molmass_ratio(tp); rtol = 1e-3)
        Test.@test 0.6 < SSCF.molmass_ratio(tp) < 0.65
    end

    # (T [K], p [Pa], q_tot [kg/kg]) across surface → upper troposphere
    pts = [(288.0, 1.0e5, 0.012), (280.0, 9.0e4, 0.008), (270.0, 7.0e4, 0.004),
           (250.0, 5.0e4, 0.001), (230.0, 3.0e4, 2.0e-4), (290.0, 1.0e5, 0.02)]

    Test.@testset "seam functions agree over the sounding" begin
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
