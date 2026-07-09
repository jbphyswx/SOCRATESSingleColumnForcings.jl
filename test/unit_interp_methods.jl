using Test: Test
using NCDatasets: NCDatasets as NC
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# Unit tests for _materialize, interp_along_dim, and var_to_new_coord.
#
# These tests use only plain Julia arrays; no NC file I/O is needed.
# They specifically cover the class of failures where NCDatasets lazy
# DiskArrays (BroadcastDiskArray with Union{Missing,Float64} elements)
# escaped materialization and caused mapslices/permutedims to crash.
# ---------------------------------------------------------------------------

Test.@testset "Interpolation method unit tests" begin

    # ------------------------------------------------------------------
    Test.@testset "_materialize: plain Array passes through unchanged" begin
        a = [1.0, 2.0, 3.0]
        Test.@test SSCF._materialize(a) === a          # exact same object (no copy)

        m = [1.0 2.0; 3.0 4.0]
        Test.@test SSCF._materialize(m) === m
    end

    # ------------------------------------------------------------------
    Test.@testset "_materialize: Union{Missing,Float64} array → Float64, missing→NaN" begin
        a_miss = Union{Missing, Float64}[1.0, missing, 3.0]
        result = SSCF._materialize(a_miss)
        Test.@test result isa Array{Float64}
        Test.@test result[1] == 1.0
        Test.@test isnan(result[2])
        Test.@test result[3] == 3.0
    end

    # ------------------------------------------------------------------
    Test.@testset "_materialize: 2D Union{Missing,Float64} matrix" begin
        m = Union{Missing, Float64}[1.0 missing; 3.0 4.0]
        rm = SSCF._materialize(m)
        Test.@test rm isa Matrix{Float64}
        Test.@test rm[1, 1] == 1.0
        Test.@test isnan(rm[1, 2])
        Test.@test rm[2, 1] == 3.0
        Test.@test rm[2, 2] == 4.0
    end

    # ------------------------------------------------------------------
    Test.@testset "_materialize: SubArray becomes a plain Array" begin
        parent = [1.0, 2.0, 3.0, 4.0]
        sub = view(parent, 2:4)
        result = SSCF._materialize(sub)
        Test.@test result isa Array
        Test.@test result == [2.0, 3.0, 4.0]
    end

    # ------------------------------------------------------------------
    Test.@testset "_materialize: lazy broadcast preserving element type" begin
        # Simulate the pattern from forcings.jl:
        #   subsidence = -(ω .- dpdt .* f_p) ./ (ρ .* g)
        # where ω is a Union{Missing,Float64} array.
        omega   = Union{Missing, Float64}[1.0, 2.0, missing, 4.0]
        rho     = [1.0, 1.0, 1.0, 1.0]
        lazy    = omega ./ rho          # a Broadcasted / lazy result
        result  = SSCF._materialize(lazy)
        Test.@test result isa Array{Float64}
        Test.@test result[1] == 1.0
        Test.@test result[2] == 2.0
        Test.@test isnan(result[3])
        Test.@test result[4] == 4.0
    end

    # ------------------------------------------------------------------
    #   interp_along_dim: basic correctness with plain Float64 array
    # ------------------------------------------------------------------
    Test.@testset "interp_along_dim: 2D (z × t) linear interpolation, full-array coord" begin
        # var[z, t]: 3 levels × 2 time steps; values increase linearly with z
        var  = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]
        z_in = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]   # same shape as var
        z_out = [1.5, 2.5]

        result = SSCF.interp_along_dim(
            var, 1, z_in;
            interp_dim_out = z_out,
            interp_dim_in_is_full_array = true,
        )
        Test.@test size(result) == (2, 2)
        Test.@test isapprox(result[1, 1], 15.0; atol = 1e-10)   # t=1, z=1.5
        Test.@test isapprox(result[1, 2], 45.0; atol = 1e-10)   # t=2, z=1.5
        Test.@test isapprox(result[2, 1], 25.0; atol = 1e-10)   # t=1, z=2.5
        Test.@test isapprox(result[2, 2], 55.0; atol = 1e-10)   # t=2, z=2.5
    end

    Test.@testset "interp_along_dim: 2D (z × t) linear interpolation, 1D coord" begin
        var   = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]
        z_in  = [1.0, 2.0, 3.0]   # 1D — not same shape → interp_dim_in_is_full_array = false
        z_out = [1.5, 2.5]

        result = SSCF.interp_along_dim(
            var, 1, z_in;
            interp_dim_out = z_out,
            interp_dim_in_is_full_array = false,
        )
        Test.@test size(result) == (2, 2)
        Test.@test isapprox(result[1, 1], 15.0; atol = 1e-10)
        Test.@test isapprox(result[2, 2], 55.0; atol = 1e-10)
    end

    # ------------------------------------------------------------------
    #   interp_along_dim: Missing-typed input must not crash
    # ------------------------------------------------------------------
    Test.@testset "interp_along_dim: Union{Missing,Float64} input does not crash" begin
        var  = Union{Missing, Float64}[10.0 40.0; missing 50.0; 30.0 60.0]
        z_in = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out = [1.5, 2.5]

        Test.@test_nowarn SSCF.interp_along_dim(
            var, 1, z_in;
            interp_dim_out = z_out,
            interp_dim_in_is_full_array = true,
        )
        result = SSCF.interp_along_dim(
            var, 1, z_in;
            interp_dim_out = z_out,
            interp_dim_in_is_full_array = true,
        )
        Test.@test eltype(result) <: AbstractFloat
        Test.@test !any(ismissing, result)
    end

    # ------------------------------------------------------------------
    #   var_to_new_coord: unweighted round-trip
    # ------------------------------------------------------------------
    Test.@testset "var_to_new_coord: unweighted basic interpolation" begin
        var   = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]
        z_in  = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out = [1.5, 2.5]

        result = SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out)
        Test.@test size(result) == (2, 2)
        Test.@test isapprox(result[1, 1], 15.0; atol = 1e-10)
        Test.@test isapprox(result[2, 2], 55.0; atol = 1e-10)
    end

    # ------------------------------------------------------------------
    #   var_to_new_coord: density-weighted gives identical result for constant weight
    # ------------------------------------------------------------------
    Test.@testset "var_to_new_coord: uniform weight equals unweighted" begin
        var    = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]
        weight = Float64[1.0 1.0; 1.0 1.0; 1.0 1.0]   # uniform → result must equal unweighted
        z_in   = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out  = [1.5, 2.5]

        result_w  = SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out, weight)
        result_uw = SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out)
        Test.@test isapprox(result_w, result_uw; atol = 1e-10)
    end

    # ------------------------------------------------------------------
    #   var_to_new_coord: Union{Missing,Float64} var + weight must not crash
    #   (this is the direct reproduction of the subsidence BroadcastDiskArray crash)
    # ------------------------------------------------------------------
    Test.@testset "var_to_new_coord: Union{Missing,Float64} weighted does not crash" begin
        var    = Union{Missing, Float64}[10.0 40.0; 20.0 50.0; 30.0 60.0]
        weight = Union{Missing, Float64}[1.0 1.0; 1.0 1.0; 1.0 1.0]
        z_in   = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out  = [1.5, 2.5]

        Test.@test_nowarn SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out, weight)
        result = SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out, weight)
        Test.@test eltype(result) <: AbstractFloat
        Test.@test !any(ismissing, result)
        Test.@test isapprox(result[1, 1], 15.0; atol = 1e-10)
    end

    # ------------------------------------------------------------------
    #   var_to_new_coord: lazy-broadcast input (simulates DiskArray arithmetic)
    # ------------------------------------------------------------------
    Test.@testset "var_to_new_coord: lazy broadcast (simulated BroadcastDiskArray) does not crash" begin
        # Simulate the production code pattern where omega has one flagged (SSCF.Interpolation.FastLinear1DInterpolationissing) measurement:
        #   measurment at z=2, t=1 is missing  →  arithmetic keeps Union{Missing,Float64} type
        omega_nc = Union{Missing, Float64}[10.0 40.0; missing 50.0; 30.0 60.0]
        rho_nc   = Union{Missing, Float64}[1.2 1.2; 1.1 1.1; 1.0 1.0]
        # -omega / (rho * g): one missing element keeps the element type as Union{Missing,Float64}
        lazy_sub = -(omega_nc ./ (rho_nc .* 9.81))
        Test.@test eltype(lazy_sub) >: Missing        # confirms Union{Missing,Float64}
        z_in = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out = [1.5, 2.5]

        Test.@test_nowarn SSCF.var_to_new_coord(lazy_sub, z_in, 1; coord_new = z_out)
        result = SSCF.var_to_new_coord(lazy_sub, z_in, 1; coord_new = z_out)
        Test.@test eltype(result) <: AbstractFloat
        Test.@test !any(ismissing, result)
    end

end

# ---------------------------------------------------------------------------
# Unit tests for the interpolation primitives in the `Interpolation` submodule:
# `build_spline` (backing-generic + eltype promotion + mixed backing),
# `drop_collinear_nodes` (exactness + Vector/SVector backing), and `interpolate_1d`.
# ---------------------------------------------------------------------------
Test.@testset "Interpolation primitives" begin
    ext = SSCF.Interpolation.ExtrapolateBoundaryCondition()
    Test.@testset "build_spline: Vector backing, linear correctness" begin
        s = SSCF.Interpolation.build_spline(SSCF.Interpolation.FastLinear1DInterpolation, [0.0, 1.0, 2.0], [0.0, 10.0, 20.0]; bc = ext)
        Test.@test s.xp isa Vector{Float64}
        Test.@test isapprox(s(0.5), 5.0; atol = 1e-12)
        Test.@test isapprox(s(1.5), 15.0; atol = 1e-12)
        Test.@test isapprox.(s.([0.5, 1.5]), [5.0, 15.0]) |> all   # broadcastable
    end

    Test.@testset "build_spline: mixed backing (Vector coord + lazy ReshapedArray data)" begin
        # fp is a non-contiguous view (ReshapedArray), a different backing than the Vector coord
        A = [0.0 99.0; 10.0 99.0; 20.0 99.0]     # (3,2); column 1 is the data
        fp_lazy = vec(view(A, :, 1))              # SubArray/ReshapedArray, ≠ Vector{Float64}
        xp = [0.0, 1.0, 2.0]
        Test.@test typeof(fp_lazy) !== typeof(xp) # genuinely different backing
        Test.@test_nowarn SSCF.Interpolation.build_spline(SSCF.Interpolation.FastLinear1DInterpolation, xp, fp_lazy; bc = ext)
        s = SSCF.Interpolation.build_spline(SSCF.Interpolation.FastLinear1DInterpolation, xp, fp_lazy; bc = ext)
        Test.@test isapprox(s(0.5), 5.0; atol = 1e-12)
        Test.@test isapprox(s(2.5), 25.0; atol = 1e-12)   # extrapolate
    end

    Test.@testset "build_spline: mixed node/value eltypes evaluate correctly (Int fp, Float64 nodes)" begin
        # nodes and values keep independent backings/eltypes; evaluation promotes, so the result is a
        # Float even when fp is stored as Int (no build-time eltype coercion needed)
        s = SSCF.Interpolation.build_spline(SSCF.Interpolation.FastLinear1DInterpolation, [0.0, 1.0, 2.0], [0, 10, 20]; bc = ext)
        Test.@test eltype(s.xp) === Float64
        Test.@test s(0.5) isa Float64
        Test.@test isapprox(s(0.5), 5.0; atol = 1e-12)
    end

    Test.@testset "build_spline: SVector backing preserved (no-drop)" begin
        xs = SSCF.Interpolation.create_svector([0.0, 1.0, 2.0])
        ys = SSCF.Interpolation.create_svector([0.0, 10.0, 20.0])
        s = SSCF.Interpolation.build_spline(SSCF.Interpolation.FastLinear1DInterpolation, xs, ys; bc = ext, drop_collinear = Val(false))
        Test.@test s.xp isa SSCF.Interpolation.StaticArrays.SVector
        Test.@test length(s.xp) == 3
        Test.@test isapprox(s(1.5), 15.0; atol = 1e-12)
    end

    Test.@testset "drop_collinear_nodes: exact prune (Vector)" begin
        # first 3 points collinear (slope 10), then a bend
        xp = [0.0, 1.0, 2.0, 3.0, 4.0]
        fp = [0.0, 10.0, 20.0, 25.0, 30.0]
        xk, yk = SSCF.Interpolation.drop_collinear_nodes(xp, fp)
        # x=1 collinear with (0,2) and x=3 collinear with (2,4) → both dropped; only the bend at x=2 kept
        Test.@test xk == [0.0, 2.0, 4.0]
        Test.@test yk == [0.0, 20.0, 30.0]
        # pruned spline evaluates identically to the unpruned one
        s_drop = SSCF.Interpolation.build_spline(SSCF.Interpolation.FastLinear1DInterpolation, xp, fp; bc = ext, drop_collinear = Val(true))
        s_full = SSCF.Interpolation.build_spline(SSCF.Interpolation.FastLinear1DInterpolation, xp, fp; bc = ext, drop_collinear = Val(false))
        Test.@test length(s_drop.xp) < length(s_full.xp)
        Test.@test all(isapprox(s_drop(x), s_full(x); atol = 1e-12) for x in 0.0:0.1:4.0)
    end

    Test.@testset "drop_collinear_nodes: all-collinear keeps only endpoints" begin
        xk, yk = SSCF.Interpolation.drop_collinear_nodes([0.0, 1.0, 2.0, 3.0], [0.0, 10.0, 20.0, 30.0])
        Test.@test xk == [0.0, 3.0]
        Test.@test yk == [0.0, 30.0]
    end

    Test.@testset "drop_collinear_nodes: SVector prunes to SVector{keep}" begin
        xs = SSCF.Interpolation.create_svector([0.0, 1.0, 2.0, 3.0])   # all collinear with ys below
        ys = SSCF.Interpolation.create_svector([0.0, 10.0, 20.0, 30.0])
        xk, yk = SSCF.Interpolation.drop_collinear_nodes(xs, ys)
        Test.@test xk isa SSCF.Interpolation.StaticArrays.SVector
        Test.@test yk isa SSCF.Interpolation.StaticArrays.SVector
        Test.@test length(xk) == 2 && xk[1] == 0.0 && xk[2] == 3.0
    end

    # Regression: each side is pruned in its OWN backing. A mismatched pair used to fall through to the
    # `similar` path and demote BOTH to an allocating `Vector`; the value backing must now be preserved.
    Test.@testset "drop_collinear_nodes: mixed (Vector coord, SVector value) preserves each backing" begin
        xp = [0.0, 1.0, 2.0, 3.0, 4.0]                                      # x=1,3 collinear → dropped
        fp = SSCF.Interpolation.create_svector([0.0, 10.0, 20.0, 25.0, 30.0])
        xk, yk = SSCF.Interpolation.drop_collinear_nodes(xp, fp)
        Test.@test xk isa Vector                                            # coord backing kept
        Test.@test yk isa SSCF.Interpolation.StaticArrays.SVector           # value backing kept (NOT demoted)
        Test.@test length(xk) == 3 && length(yk) == 3
        Test.@test xk == [0.0, 2.0, 4.0]
        Test.@test collect(yk) == [0.0, 20.0, 30.0]
    end

    Test.@testset "drop_collinear_nodes: mixed (SVector coord, Vector value) preserves each backing" begin
        xs = SSCF.Interpolation.create_svector([0.0, 1.0, 2.0, 3.0, 4.0])
        fv = [0.0, 10.0, 20.0, 25.0, 30.0]
        xk, yk = SSCF.Interpolation.drop_collinear_nodes(xs, fv)
        Test.@test xk isa SSCF.Interpolation.StaticArrays.SVector
        Test.@test yk isa Vector
        Test.@test length(xk) == 3 && length(yk) == 3
    end

    Test.@testset "drop_collinear_nodes: range coord stays a range when pruned nodes are uniform" begin
        xr = 0.0:1.0:3.0                       # uniform grid
        fp = [0.0, 10.0, 20.0, 30.0]           # all collinear → only endpoints kept (2 nodes, trivially uniform)
        xk, yk = SSCF.Interpolation.drop_collinear_nodes(xr, fp)
        Test.@test xk isa AbstractRange        # range backing preserved (no Vector demotion)
        Test.@test collect(xk) == [0.0, 3.0]
    end

    Test.@testset "drop_collinear_nodes: range coord errors when pruned nodes are non-uniform" begin
        xr = 0.0:1.0:4.0
        # kept indices 1,2,5 → pruned x = [0,1,4] (gaps 1 then 3) → not representable as a range
        fp = [5.0, 10.0, 20.0, 30.0, 40.0]
        Test.@test_throws ErrorException SSCF.Interpolation.drop_collinear_nodes(xr, fp)
    end

    Test.@testset "interpolate_1d: one-shot build + evaluate" begin
        y = SSCF.Interpolation.interpolate_1d([0.5, 1.5], [0.0, 1.0, 2.0], [0.0, 10.0, 20.0], SSCF.Interpolation.FastLinear1DInterpolation; bc = ext)
        Test.@test isapprox(y, [5.0, 15.0]; atol = 1e-12)
    end

    Test.@testset "safe_integrate: trapezoidal-exact for a piecewise-linear spline over arbitrary ranges" begin
        # ∫ of a piecewise-linear function equals the trapezoid rule exactly, for limits on or between nodes
        z = collect(0.0:100.0:1000.0)
        f = 1.0 .+ 0.001 .* z                       # exact line
        spl = SSCF.Interpolation.build_spline(SSCF.Interpolation.FastLinear1DInterpolation, z, f; bc = ext)
        line_int(a, b) = ((1 + 0.001a) + (1 + 0.001b)) / 2 * (b - a)   # trapezoid == exact for the line
        Test.@test isapprox(SSCF.Interpolation.safe_integrate(spl, 100.0, 300.0; bc = ext), line_int(100, 300); rtol = 1e-12)  # node-aligned
        Test.@test isapprox(SSCF.Interpolation.safe_integrate(spl, 125.0, 375.0; bc = ext), line_int(125, 375); rtol = 1e-12)  # lower limit mid-segment
        Test.@test isapprox(SSCF.Interpolation.safe_integrate(spl, 137.0, 862.0; bc = ext), line_int(137, 862); rtol = 1e-12)  # both limits mid-segment
        Test.@test isapprox(SSCF.Interpolation.safe_integrate(spl, 120.0, 180.0; bc = ext), line_int(120, 180); rtol = 1e-12)  # within a single segment
    end
end
