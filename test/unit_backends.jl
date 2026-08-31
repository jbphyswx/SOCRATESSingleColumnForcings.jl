using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Thermodynamics: Thermodynamics as TD
using Dierckx: Dierckx
using PCHIPInterpolation: PCHIPInterpolation
using ForwardDiff: ForwardDiff
using Interpolations: Interpolations

# ---------------------------------------------------------------------------
# Backend conformance. Every shipped backend is exercised through the SAME signatures the pipeline
# uses, so a backend that is missing one fails here rather than at a user's call site. This is the
# gate that a per-backend smoke test cannot be: `test/extensions.jl` builds each interpolation backend
# without `drop_collinear`, which is a signature the pipeline never calls, so a backend lacking that
# keyword passed there while erroring inside `interp_along_dim`.
#
# Fully qualified paths throughout, no local aliases: these files share one `Main`, and `I` collides
# with `LinearAlgebra.I`.
# ---------------------------------------------------------------------------

backend_ext(name) = Base.get_extension(SSCF, name)

# (label, method) for every interpolation backend the package ships.
function shipped_interpolation_methods()
    d = backend_ext(:SOCRATESSingleColumnForcingsDierckxExt)
    p = backend_ext(:SOCRATESSingleColumnForcingsPCHIPInterpolationExt)
    i = backend_ext(:SOCRATESSingleColumnForcingsInterpolationsExt)
    return (
        ("FastLinear", SSCF.Interpolation.FastLinear1DInterpolation),
        ("Dierckx(k=1)", d.DierckxSpline1DInterpolationMethod(1)),
        ("Dierckx(k=3)", d.DierckxSpline1DInterpolationMethod(3)),
        ("PCHIP", p.PCHIPInterpolationMethod()),
        (
            "Interpolations(Gridded(Linear))",
            i.InterpolationsInterpolationMethod(Interpolations.Gridded(Interpolations.Linear())),
        ),
    )
end

Test.@testset "Interpolation backend conformance" begin
    extrap = SSCF.Interpolation.ExtrapolateBoundaryCondition()
    xp = collect(0.0:1.0:10.0)
    fp = @. 1.0 + 0.5 * xp           # a line: every method must reproduce it exactly
    at(x) = 1.0 + 0.5 * x

    Test.@testset "$label" for (label, method) in shipped_interpolation_methods()
        # The pipeline always passes `drop_collinear`; a backend missing it errors in interp_along_dim.
        Test.@testset "build_spline accepts drop_collinear (the pipeline's signature)" begin
            for dc in (Val(false), Val(true))
                spl = SSCF.Interpolation.build_spline(method, xp, fp; bc = extrap, drop_collinear = dc)
                Test.@test isapprox(spl(2.5), at(2.5); atol = 1e-8)
                Test.@test isapprox(spl(7.25), at(7.25); atol = 1e-8)
            end
        end

        Test.@testset "evaluates in and out of range" begin
            spl = SSCF.Interpolation.build_spline(method, xp, fp; bc = extrap)
            Test.@test isapprox(spl(0.0), at(0.0); atol = 1e-8)
            Test.@test isapprox(spl(10.0), at(10.0); atol = 1e-8)
            Test.@test isfinite(spl(12.0))   # extrapolate bc: a value, not a throw
        end
    end
end

Test.@testset "safe_integrate over a line is exact where implemented" begin
    extrap = SSCF.Interpolation.ExtrapolateBoundaryCondition()
    xp = collect(0.0:1.0:10.0)
    fp = @. 1.0 + 0.5 * xp
    exact(a, b) = (b - a) + 0.25 * (b^2 - a^2)   # ∫(1 + x/2)dx

    Test.@testset "$label" for (label, method) in shipped_interpolation_methods()
        spl = SSCF.Interpolation.build_spline(method, xp, fp; bc = extrap)
        if hasmethod(SSCF.Interpolation.safe_integrate, Tuple{typeof(spl), Float64, Float64})
            Test.@test isapprox(
                SSCF.Interpolation.safe_integrate(spl, 1.0, 3.0; bc = extrap), exact(1.0, 3.0); rtol = 1e-8,
            )
            Test.@test isapprox(
                SSCF.Interpolation.safe_integrate(spl, 0.0, 10.0; bc = extrap), exact(0.0, 10.0); rtol = 1e-8,
            )
            Test.@test isapprox(
                SSCF.Interpolation.safe_integrate(spl, 3.0, 1.0; bc = extrap),
                -SSCF.Interpolation.safe_integrate(spl, 1.0, 3.0; bc = extrap);
                rtol = 1e-8,
            )
        else
            # Recorded, not silently skipped: conservative regridding cannot work for this backend.
            @info "safe_integrate not implemented for $label; conservative regridding is unavailable there"
            Test.@test_broken hasmethod(
                SSCF.Interpolation.safe_integrate, Tuple{typeof(spl), Float64, Float64},
            )
        end
    end
end

Test.@testset "create_bc reports invalid names" begin
    for good in (:extrapolate, :error, :nearest)
        Test.@test SSCF.Interpolation.create_bc(good) isa SSCF.Interpolation.ValidBoundaryConditions
        Test.@test SSCF.Interpolation.create_bc(String(good)) isa SSCF.Interpolation.ValidBoundaryConditions
    end
    # The Symbol/String funnel must reach the helpful message, not a bare MethodError.
    err = try
        SSCF.Interpolation.create_bc(:cubic)
        nothing
    catch e
        e
    end
    Test.@test err isa ErrorException
    Test.@test occursin("cubic", sprint(showerror, err))
end

# ---------------------------------------------------------------------------
# Thermodynamics backends: both must answer every generic function the pipeline calls on them.
# ---------------------------------------------------------------------------

function accurate_thermo_params()
    b = SSCF.DefaultThermodynamicsBackend()
    R_d = SSCF.R_d(b)
    R_v = SSCF.R_v(b)
    cp_d = SSCF.cp_d(b)
    Ru = 8.3144598  # [J/mol/K]; molmass_{dryair,water} = Ru/R_{d,v} reproduce R_d, R_v
    vals = (;
        T_0 = 273.16, T_triple = 273.16, T_freeze = SSCF.T_freeze(b), T_icenuc = SSCF.T_icenuc(b),
        T_min = 1.0, T_max = 1000.0, T_init_min = 150.0, T_surf_ref = 290.0, T_min_ref = 220.0,
        entropy_reference_temperature = 298.15, MSLP = 101325.0, p_ref_theta = SSCF.p_ref(b),
        press_triple = 611.657, R_d = R_d, R_v = R_v, cp_d = cp_d, cp_v = SSCF.cp_v(b),
        cp_l = SSCF.cp_l(b), cp_i = SSCF.cp_i(b), LH_v0 = SSCF.L_v0(b), LH_s0 = SSCF.L_s0(b),
        entropy_dry_air = 6864.8, entropy_water_vapor = 10513.6, grav = SSCF.grav(b), pow_icenuc = 1.0,
        q_min = eps(Float64),
        gas_constant = Ru, molmass_dryair = Ru / R_d, molmass_water = Ru / R_v, kappa_d = R_d / cp_d,
    )
    fns = fieldnames(TD.Parameters.ThermodynamicsParameters)
    absent = filter(fn -> !haskey(vals, fn), fns)
    # Fail loudly rather than filling unknown fields with NaN: a field the pipeline reads would
    # otherwise silently produce NaN results.
    isempty(absent) ||
        error("ThermodynamicsParameters gained fields this fixture does not set: $(absent)")
    return TD.Parameters.ThermodynamicsParameters{Float64}(; (fn => vals[fn] for fn in fns)...)
end

Test.@testset "Thermodynamics backend conformance" begin
    naive = SSCF.DefaultThermodynamicsBackend()
    accurate = accurate_thermo_params()
    T, p, q = 275.0, 9.0e4, 4.0e-3
    pg, Tg = 1.0e5, 288.0

    Test.@testset "$label answers every pipeline call" for (label, b) in
                                                          (("default", naive), ("Thermodynamics.jl", accurate))
        Test.@test SSCF.R_d(b) > 0
        Test.@test SSCF.R_v(b) > SSCF.R_d(b)
        Test.@test SSCF.grav(b) > 0
        Test.@test 0 < SSCF.molmass_ratio(b) < 1     # ε = M_v/M_d ≈ 0.622, NOT Thermodynamics' inverse
        Test.@test SSCF.R_d(b, Float32) isa Float32  # the 2-arg form the pipeline uses to fix precision

        Test.@test SSCF.dry_pottemp(b, T, p) > T      # p < p_ref
        Test.@test SSCF.virtual_temperature(b, T, p, q) >= T
        Test.@test SSCF.air_density(b, T, p, q) > 0
        Test.@test SSCF.air_density(b, T, p, q, 0.0, 0.0) > 0
        Test.@test SSCF.liquid_ice_pottemp(b, T, p, q) > 0
        Test.@test SSCF.liquid_ice_pottemp(b, T, p, q, 0.0, 0.0) > 0
        Test.@test 0 < SSCF.q_vap_saturation_liq(b, T, p) < 1
        Test.@test 0 < SSCF.saturation_specific_humidity_from_pT(b, pg, Tg) < 1
        Test.@test 0 < SSCF.saturation_mixing_ratio_from_pT(b, pg, Tg) < 1

        c = SSCF.equilibrium_condensate(b, T, p, q)
        Test.@test c.q_liq >= 0 && c.q_ice >= 0
        Test.@test c.q_liq + c.q_ice <= q

        s = SSCF.saturation_adjust_pθq(b, p, SSCF.dry_pottemp(b, T, p), q)
        Test.@test isapprox(s.T, T; atol = 0.2)       # unsaturated round-trip recovers T
        Test.@test s.q_liq >= 0 && s.q_ice >= 0
    end

    Test.@testset "enforced liquid fraction" begin
        # λ is the whole point of the knob: forcing all-liquid must suppress ice below freezing.
        Tc, pc, qc = 260.0, 8.0e4, 8.0e-3            # supercooled and saturated enough to condense
        ramp = SSCF.equilibrium_condensate(naive, Tc, pc, qc)
        allliq = SSCF.equilibrium_condensate(naive, Tc, pc, qc; λ = 1.0)
        Test.@test ramp.q_ice > 0                     # the default ramp puts some condensate in ice here
        Test.@test allliq.q_ice == 0
        Test.@test allliq.q_liq > 0

        # and it must reach through the saturation adjustment, which is where the pipeline uses it
        θ = SSCF.liquid_ice_pottemp(naive, Tc, pc, qc)
        Test.@test SSCF.saturation_adjust_pθq(naive, pc, θ, qc; λ = 1.0).q_ice == 0

        # the default sentinel must NOT poison the result
        Test.@test all(isfinite, values(SSCF.saturation_adjust_pθq(naive, pc, θ, qc)))

        Test.@test SSCF.saturation_adjust_pθq(accurate, pc, θ, qc; λ = 1.0).q_ice == 0
        # both backends agree on the all-liquid partition the pipeline asks for
        Test.@test isapprox(
            SSCF.saturation_adjust_pθq(naive, pc, θ, qc; λ = 1.0).T,
            SSCF.saturation_adjust_pθq(accurate, pc, θ, qc; λ = 1.0).T;
            rtol = 1e-2,
        )
    end
end
