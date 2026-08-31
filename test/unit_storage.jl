using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# Interpolant node/value storage. Every backing the README and `docs/src/interpolation.md` advertise
# is built here through `coerce_vector` — the same entry point `get_column_forcing` uses — and then
# evaluated, so a backing that cannot represent this package's own axes fails here.
#
# ---------------------------------------------------------------------------

Test.@testset "Storage backings" begin
    t_int = collect(0:3600:(3600 * 10))          # the Atlas hourly axis: integer seconds, exactly uniform
    t_float = Float64.(t_int)
    f = @. 2.0 + 3.0e-4 * t_float

    Test.@testset "UniformRange on an integer axis" begin
        # `inv(step)` is not an integer, so a UniformRange whose inv_step shares the node element
        # type cannot exist for this axis at all.
        r = SSCF.Interpolation.coerce_vector(Tuple{SSCF.Interpolation.UniformRange, Nothing}, t_int)
        Test.@test r isa SSCF.Interpolation.UniformRange
        Test.@test length(r) == length(t_int)
        Test.@test first(r) == first(t_int)
        Test.@test last(r) == last(t_int)
        Test.@test step(r) == 3600
        Test.@test all(r[i] == t_int[i] for i in eachindex(t_int))

        # and it must evaluate like the plain-Vector interpolant it is an optimization of
        # (a range coordinate requires drop_collinear = Val(false); pruning cannot preserve a range)
        itp_r = SSCF.Interpolation.build_spline(
            SSCF.Interpolation.FastLinear1DInterpolation, r, f;
            bc = SSCF.Interpolation.ExtrapolateBoundaryCondition(), drop_collinear = Val(false),
        )
        itp_v = SSCF.Interpolation.build_spline(
            SSCF.Interpolation.FastLinear1DInterpolation, t_float, f;
            bc = SSCF.Interpolation.ExtrapolateBoundaryCondition(),
        )
        for x in (0.0, 1234.5, 3600.0, 18_000.0, 36_000.0)
            Test.@test isapprox(itp_r(x), itp_v(x); rtol = 1e-12)
        end
    end

    Test.@testset "UniformRange on a float axis, ascending and descending" begin
        up = SSCF.Interpolation.coerce_vector(
            Tuple{SSCF.Interpolation.UniformRange, Nothing}, collect(0.0:0.5:5.0),
        )
        Test.@test step(up) == 0.5 && last(up) == 5.0
        down = SSCF.Interpolation.coerce_vector(
            Tuple{SSCF.Interpolation.UniformRange, Nothing}, collect(5.0:-0.5:0.0),
        )
        Test.@test step(down) == -0.5 && last(down) == 0.0
    end

    Test.@testset "UniformRange rejects a non-uniform axis" begin
        Test.@test_throws ErrorException SSCF.Interpolation.coerce_vector(
            Tuple{SSCF.Interpolation.UniformRange, Nothing}, [0.0, 1.0, 3.0],
        )
    end

    Test.@testset "each advertised backing round-trips through coerce_vector" begin
        Test.@test SSCF.Interpolation.coerce_vector(Tuple{Vector, Float32}, t_float) isa Vector{Float32}
        Test.@test SSCF.Interpolation.coerce_vector(Tuple{Vector, Nothing}, t_int) isa Vector{Int}
        Test.@test SSCF.Interpolation.coerce_vector(Tuple{Nothing, Float32}, t_float) isa AbstractVector{Float32}
        Test.@test SSCF.Interpolation.coerce_vector(Tuple{Nothing, Nothing}, t_float) === t_float # full passthrough
        Test.@test SSCF.Interpolation.coerce_vector(Tuple{SSCF.StaticArrays.SVector, Float64}, t_float) isa
                   SSCF.StaticArrays.SVector
        Test.@test SSCF.Interpolation.coerce_vector(Tuple{StepRangeLen, Nothing}, t_int) isa AbstractRange
    end

    Test.@testset "Constant / ConstantVector need all-equal values" begin
        flat = fill(7.5, 6)
        c = SSCF.Interpolation.coerce_vector(Tuple{SSCF.Interpolation.Constant, Nothing}, flat)
        Test.@test c isa SSCF.Interpolation.Constant
        Test.@test length(c) == 1 && c[1] == 7.5

        cv = SSCF.Interpolation.coerce_vector(Tuple{SSCF.Interpolation.ConstantVector, Nothing}, flat)
        Test.@test cv isa SSCF.Interpolation.ConstantVector
        Test.@test length(cv) == length(flat)
        Test.@test all(cv[i] == 7.5 for i in eachindex(flat))

        Test.@test_throws ErrorException SSCF.Interpolation.coerce_vector(
            Tuple{SSCF.Interpolation.Constant, Nothing}, [1.0, 2.0],
        )
        Test.@test_throws ErrorException SSCF.Interpolation.coerce_vector(
            Tuple{SSCF.Interpolation.ConstantVector, Nothing}, [1.0, 2.0],
        )
    end

    Test.@testset "a constant-valued interpolant evaluates to its value" begin
        xs = collect(0.0:1.0:5.0)
        itp = SSCF.Interpolation.build_spline(
            SSCF.Interpolation.FastLinear1DInterpolation,
            xs,
            SSCF.Interpolation.coerce_vector(
                Tuple{SSCF.Interpolation.ConstantVector, Nothing}, fill(4.25, length(xs)),
            );
            bc = SSCF.Interpolation.ExtrapolateBoundaryCondition(),
        )
        Test.@test itp(0.0) == 4.25
        Test.@test itp(2.5) == 4.25
        Test.@test itp(9.0) == 4.25   # extrapolate bc: still the constant
    end
end
