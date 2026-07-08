using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Thermodynamics: Thermodynamics
using Dierckx: Dierckx
using PCHIPInterpolation: PCHIPInterpolation
using ForwardDiff: ForwardDiff
using Interpolations: Interpolations
using NonNegLeastSquares: NonNegLeastSquares

# ---------------------------------------------------------------------------
# Every weakdep extension must precompile and activate once its trigger package is loaded, and the
# interpolation backend it provides must build & evaluate. (Extensions that are never loaded never
# precompile, so a load-time bug in one stays invisible until it is actually triggered.)
# ---------------------------------------------------------------------------
Test.@testset "Extensions" begin
    I = SSCF.Interpolation
    extrap = I.ExtrapolateBoundaryCondition()
    xp = collect(0.0:1.0:10.0)
    fp = @. 1.0 + 0.5 * xp                    # a line: every method must reproduce it exactly
    at(x) = 1.0 + 0.5 * x

    Test.@testset "all weakdep extensions active" begin
        for name in (
            :SOCRATESSingleColumnForcingsThermodynamicsExt,
            :SOCRATESSingleColumnForcingsDierckxExt,
            :SOCRATESSingleColumnForcingsPCHIPInterpolationExt,
            :SOCRATESSingleColumnForcingsInterpolationsExt,
            :SOCRATESSingleColumnForcingsNonNegLeastSquaresExt,
        )
            Test.@test Base.get_extension(SSCF, name) !== nothing
        end
    end

    Test.@testset "Dierckx backend builds & evaluates" begin
        dext = Base.get_extension(SSCF, :SOCRATESSingleColumnForcingsDierckxExt)
        spl = I.build_spline(dext.DierckxSpline1DInterpolationMethod(3), xp, fp; bc = extrap)
        Test.@test isapprox(spl(2.5), at(2.5); atol = 1e-8)
    end

    Test.@testset "PCHIP backend builds & evaluates" begin
        pext = Base.get_extension(SSCF, :SOCRATESSingleColumnForcingsPCHIPInterpolationExt)
        spl = I.build_spline(pext.PCHIPInterpolationMethod(), xp, fp; bc = extrap)
        Test.@test isapprox(spl(2.5), at(2.5); atol = 1e-8)
    end

    Test.@testset "Interpolations backend builds & evaluates" begin
        iext = Base.get_extension(SSCF, :SOCRATESSingleColumnForcingsInterpolationsExt)
        m = iext.InterpolationsInterpolationMethod(Interpolations.Gridded(Interpolations.Linear()))
        spl = I.build_spline(m, xp, fp; bc = extrap)
        Test.@test isapprox(spl(2.5), at(2.5); atol = 1e-8)
    end
end
