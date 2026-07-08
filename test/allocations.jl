using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# The interpolant evaluation is the TC per-timestep hot path and must be allocation-free for every
# backing (SVector = isbits/stack, and Vector on both uniform O(1) and irregular search paths). The
# function barrier `_eval` keeps the measurement free of global-scope boxing; each interpolant is
# warmed up (compiled) before `@allocated`.
# ---------------------------------------------------------------------------
_eval(itp, x) = itp(x)

Test.@testset "hot-path eval is allocation-free" begin
    I = SSCF.Interpolation
    ext = I.ExtrapolateBoundaryCondition()

    su = I.build_spline(I.FastLinear1DInterpolation, collect(0.0:1.0:40.0), sin.(0.0:1.0:40.0); bc = ext, drop_collinear = false)          # Vector, uniform (O(1))
    si = I.build_spline(I.FastLinear1DInterpolation, [0.0, 1.0, 3.0, 7.0, 15.0, 31.0], cos.([0.0, 1.0, 3.0, 7.0, 15.0, 31.0]); bc = ext, drop_collinear = false)  # Vector, irregular (search)
    ss = I.build_spline(I.FastLinear1DInterpolation, I.create_svector(collect(0.0:1.0:9.0)), I.create_svector(collect(0.0:2.0:18.0)); bc = ext, drop_collinear = false)  # SVector, isbits

    _eval(su, 7.3); _eval(si, 5.0); _eval(ss, 2.5)   # warm up / compile

    Test.@test @allocated(_eval(su, 7.3)) == 0
    Test.@test @allocated(_eval(si, 5.0)) == 0
    Test.@test @allocated(_eval(ss, 2.5)) == 0
end
