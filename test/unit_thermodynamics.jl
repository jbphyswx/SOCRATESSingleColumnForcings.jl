using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# Unit tests for the default `DefaultThermodynamicsBackend`:
# Each testset pins a physical property of the fallback backend.
# ---------------------------------------------------------------------------
Test.@testset "Thermodynamics (default backend)" begin
    b = SSCF.DefaultThermodynamicsBackend()

    Test.@testset "physical constants (1-arg and FT-parameterized)" begin
        Test.@test SSCF.R_d(b) ≈ 287.058
        Test.@test SSCF.R_v(b) > SSCF.R_d(b)
        Test.@test SSCF.grav(b) ≈ 9.80665
        Test.@test 0 < SSCF.molmass_ratio(b) < 1          # ε = M_v/M_d ≈ 0.622
        Test.@test SSCF.R_d(b, Float32) isa Float32
        Test.@test SSCF.R_d(b, Float32) ≈ Float32(287.058)
    end

    Test.@testset "saturation vapor pressure over liquid" begin
        T_fr = SSCF.T_freeze(b)
        Test.@test SSCF.saturation_vapor_pressure_liquid(b, T_fr) ≈ 611.2 rtol = 1e-6  # Clausius–Clapeyron anchor
        Test.@test SSCF.saturation_vapor_pressure_liquid(b, 300.0) >
                   SSCF.saturation_vapor_pressure_liquid(b, 280.0)                     # increasing in T
    end

    Test.@testset "liquid_fraction: 1 warm, 0 cold, linear between" begin
        Test.@test SSCF.liquid_fraction(b, 300.0) == 1.0
        Test.@test SSCF.liquid_fraction(b, 200.0) == 0.0
        Test.@test 0.0 < SSCF.liquid_fraction(b, 250.0) < 1.0
    end

    Test.@testset "q_vap_saturation: positive; unsaturatable at low pressure" begin
        # normal conditions: finite positive saturation specific humidity
        qsat = SSCF.q_vap_saturation(b, 290.0, 90000.0)
        Test.@test qsat > 0 && qsat < 0.1
        # when e_sat ≥ p (very low pressure) the air cannot saturate: q* ≥ 1, never negative
        Test.@test SSCF.q_vap_saturation(b, 270.0, 100.0) ≥ 1.0
        Test.@test SSCF.q_vap_saturation_liquid(b, 270.0, 100.0) ≥ 1.0
    end

    Test.@testset "equilibrium_condensate: non-negative, bounded by q_tot" begin
        # dry / unsaturated: no condensate
        c0 = SSCF.equilibrium_condensate(b, 290.0, 90000.0, 1e-4)
        Test.@test c0.q_liq == 0.0 && c0.q_ice == 0.0
        # saturated: condensate positive but ≤ q_tot
        cs = SSCF.equilibrium_condensate(b, 290.0, 90000.0, 0.05)
        Test.@test cs.q_liq > 0.0
        Test.@test cs.q_liq + cs.q_ice ≤ 0.05 + 1e-12
        # at very low pressure (unsaturatable) no condensate forms even at saturation-scale q_tot
        clp = SSCF.equilibrium_condensate(b, 270.0, 100.0, 4.06e-6)
        Test.@test clp.q_liq ≈ 0.0 atol = 1e-12
        Test.@test clp.q_ice ≈ 0.0 atol = 1e-12
    end

    Test.@testset "virtual_temperature: ≈T dry, >T moist, always positive" begin
        Test.@test SSCF.virtual_temperature(b, 290.0, 90000.0, 0.0) ≈ 290.0
        Test.@test SSCF.virtual_temperature(b, 290.0, 90000.0, 0.01) > 290.0
        Test.@test SSCF.virtual_temperature(b, 270.0, 100.0, 4.06e-6) > 0.0   # stays physical at low p
    end

    Test.@testset "air_density: ideal-gas; condensate-aware form differs from vapor-only" begin
        ρ_dry = SSCF.air_density(b, 290.0, 90000.0, 0.0)
        Test.@test isapprox(ρ_dry, 90000.0 / (SSCF.R_d(b) * 290.0); rtol = 1e-6)
        # 6-arg (explicit condensate) ≠ 3-arg (all q_tot treated as vapor) when condensate is present
        q_tot = 0.02; q_liq = 0.005
        Test.@test SSCF.air_density(b, 290.0, 90000.0, q_tot) !=
                   SSCF.air_density(b, 290.0, 90000.0, q_tot, q_liq, 0.0)
    end

    Test.@testset "saturation_adjust_pθq: recovers T from θ_liq_ice (unsaturated round-trip)" begin
        # for an unsaturated parcel θ_liq_ice = dry potential temperature, so adjusting
        # (p, dry_pottemp(T), q_tot) returns the original T with no condensate
        T, p, q_tot = 285.0, 90000.0, 1e-4
        θ = SSCF.dry_pottemp(b, T, p)
        T_rec, q_liq, q_ice = SSCF.saturation_adjust_pθq(b, p, θ, q_tot)
        Test.@test isapprox(T_rec, T; atol = 1e-2)   # to the bisection tolerance
        Test.@test q_liq == 0.0 && q_ice == 0.0
    end

    Test.@testset "surface saturation mixing ratio: w_s = ε·e_s/(p−e_s)" begin
        # saturation *mixing ratio* at the surface (later converted to specific humidity by the pipeline)
        pg, Tg = 1.0e5, 288.0
        ε = SSCF.molmass_ratio(b)
        e_s = SSCF.saturation_vapor_pressure_liquid(b, Tg)
        w_s = SSCF.saturation_q_tot_from_pgTg(b, pg, Tg)
        Test.@test isapprox(w_s, ε * e_s / (pg - e_s); rtol = 1e-12)   # exact mixing-ratio formula
        Test.@test 0.008 < w_s < 0.013                                # ~10.8 g/kg at 288 K / 1000 hPa
        Test.@test w_s > SSCF.calc_qg_from_pgTg(b, pg, Tg)            # mixing ratio > specific humidity
    end
end
