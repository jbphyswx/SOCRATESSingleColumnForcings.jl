using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# The interpolant evaluation is the TC per-timestep hot path and must be allocation-free for every
# backing (SVector = isbits/stack, and Vector on both uniform O(1) and irregular search paths). The
# function barrier `_eval` keeps the measurement free of global-scope boxing; each interpolant is
# warmed up (compiled) before `@allocated`.
# ---------------------------------------------------------------------------
#=
    see https://discourse.julialang.org/t/allocations-on-field-access-even-though-no-abstract-types-are-involved/126672/6
    On 1.12 things like
        Test.@test @allocated(ss(2.5)) == 0
    do not allocate, but the field accesses inside the interpolant do seem to allocate a 16 byte single allocation on 1.11
    This form does not allocate, and we believe this shows it is only a global context lowering problem
    (also defining the interpolants su, si, ss as cosnts would fix, let block does not either, etc...)
=#
function allocation_test(itp, x) # https://discourse.julialang.org/t/allocations-on-field-access-even-though-no-abstract-types-are-involved/126672/6
    @allocated itp(x)
end

Test.@testset "hot-path eval is allocation-free" begin
    I = SSCF.Interpolation
    ext = I.ExtrapolateBoundaryCondition()

    su = I.build_spline(I.FastLinear1DInterpolation, collect(0.0:1.0:40.0), sin.(0.0:1.0:40.0); bc = ext, drop_collinear = false)          # Vector, uniform (O(1))
    si = I.build_spline(I.FastLinear1DInterpolation, [0.0, 1.0, 3.0, 7.0, 15.0, 31.0], cos.([0.0, 1.0, 3.0, 7.0, 15.0, 31.0]); bc = ext, drop_collinear = false)  # Vector, irregular (search)
    ss = I.build_spline(I.FastLinear1DInterpolation, I.create_svector(collect(0.0:1.0:9.0)), I.create_svector(collect(0.0:2.0:18.0)); bc = ext, drop_collinear = false)  # SVector, isbits

    su(7.3); si(5.0); ss(2.5)   # warm up / compile
    allocation_test(su, 7.3); allocation_test(si, 5.0); allocation_test(ss, 2.5); 

    Test.@test allocation_test(su, 7.3) == 0
    Test.@test allocation_test(si, 5.0) == 0
    Test.@test allocation_test(ss, 2.5) == 0
end
