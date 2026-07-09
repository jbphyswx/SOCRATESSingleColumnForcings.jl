using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
import JET

# ---------------------------------------------------------------------------
# Type-stability of the interpolant evaluation hot path (the returned time-splines are evaluated once
# per model timestep). Each backing (uniform, irregular, SVector) must infer to a concrete
# return type with no runtime dispatch (JET `@report_opt` clean).
# ---------------------------------------------------------------------------
Test.@testset "hot-path inferrability" begin
    I = SSCF.Interpolation
    ext = I.ExtrapolateBoundaryCondition()

    xu = collect(0.0:1.0:40.0); fu = @. sin(xu / 5)              # uniform grid → UniformNodes
    su = I.build_spline(I.FastLinear1DInterpolation, xu, fu; bc = ext, drop_collinear = Val(false))
    xi = [0.0, 1.0, 3.0, 7.0, 15.0, 31.0]; fi = @. cos(xi / 4)   # irregular grid → IrregularNodes
    si = I.build_spline(I.FastLinear1DInterpolation, xi, fi; bc = ext, drop_collinear = Val(false))
    ss = I.build_spline(I.FastLinear1DInterpolation,           # SVector backing (isbits, alloc-free)
        I.create_svector(collect(0.0:1.0:4.0)), I.create_svector([0.0, 1.0, 4.0, 9.0, 16.0]);
        bc = ext, drop_collinear = Val(false))
    xr = 0.0:1.0:40.0                                          # range coord (StepRangeLen) → O(1); the storage the public forcing functions produce
    sr = I.build_spline(I.FastLinear1DInterpolation, xr, sin.(xr ./ 5); bc = ext, drop_collinear = Val(false))

    Test.@test (Test.@inferred su(7.3)) isa Float64
    Test.@test (Test.@inferred si(5.0)) isa Float64
    Test.@test (Test.@inferred ss(2.5)) isa Float64
    Test.@test (Test.@inferred sr(7.3)) isa Float64

    Test.@test isempty(JET.get_reports(JET.@report_opt su(7.3)))
    Test.@test isempty(JET.get_reports(JET.@report_opt si(5.0)))
    Test.@test isempty(JET.get_reports(JET.@report_opt ss(2.5)))
    Test.@test isempty(JET.get_reports(JET.@report_opt sr(7.3)))
end
