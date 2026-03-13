using Test
import NCDatasets as NC
import SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# Unit tests for _materialize, interp_along_dim, and var_to_new_coord.
#
# These tests use only plain Julia arrays; no NC file I/O is needed.
# They specifically cover the class of failures where NCDatasets lazy
# DiskArrays (BroadcastDiskArray with Union{Missing,Float64} elements)
# escaped materialization and caused mapslices/permutedims to crash.
# ---------------------------------------------------------------------------

@testset "Interpolation method unit tests" begin

    # ------------------------------------------------------------------
    @testset "_materialize: plain Array passes through unchanged" begin
        a = [1.0, 2.0, 3.0]
        @test SSCF._materialize(a) === a          # exact same object (no copy)

        m = [1.0 2.0; 3.0 4.0]
        @test SSCF._materialize(m) === m
    end

    # ------------------------------------------------------------------
    @testset "_materialize: Union{Missing,Float64} array → Float64, missing→NaN" begin
        a_miss = Union{Missing, Float64}[1.0, missing, 3.0]
        result = SSCF._materialize(a_miss)
        @test result isa Array{Float64}
        @test result[1] == 1.0
        @test isnan(result[2])
        @test result[3] == 3.0
    end

    # ------------------------------------------------------------------
    @testset "_materialize: 2D Union{Missing,Float64} matrix" begin
        m = Union{Missing, Float64}[1.0 missing; 3.0 4.0]
        rm = SSCF._materialize(m)
        @test rm isa Matrix{Float64}
        @test rm[1, 1] == 1.0
        @test isnan(rm[1, 2])
        @test rm[2, 1] == 3.0
        @test rm[2, 2] == 4.0
    end

    # ------------------------------------------------------------------
    @testset "_materialize: SubArray becomes a plain Array" begin
        parent = [1.0, 2.0, 3.0, 4.0]
        sub = view(parent, 2:4)
        result = SSCF._materialize(sub)
        @test result isa Array
        @test result == [2.0, 3.0, 4.0]
    end

    # ------------------------------------------------------------------
    @testset "_materialize: lazy broadcast preserving element type" begin
        # Simulate the pattern from process_case.jl:
        #   subsidence = -(ω .- dpdt .* f_p) ./ (ρ .* g)
        # where ω is a Union{Missing,Float64} array.
        omega   = Union{Missing, Float64}[1.0, 2.0, missing, 4.0]
        rho     = [1.0, 1.0, 1.0, 1.0]
        lazy    = omega ./ rho          # a Broadcasted / lazy result
        result  = SSCF._materialize(lazy)
        @test result isa Array{Float64}
        @test result[1] == 1.0
        @test result[2] == 2.0
        @test isnan(result[3])
        @test result[4] == 4.0
    end

    # ------------------------------------------------------------------
    #   interp_along_dim: basic correctness with plain Float64 array
    # ------------------------------------------------------------------
    @testset "interp_along_dim: 2D (z × t) linear interpolation, full-array coord" begin
        # var[z, t]: 3 levels × 2 time steps; values increase linearly with z
        var  = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]
        z_in = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]   # same shape as var
        z_out = [1.5, 2.5]

        result = SSCF.interp_along_dim(
            var, 1, z_in;
            interp_dim_out = z_out,
            interp_dim_in_is_full_array = true,
        )
        @test size(result) == (2, 2)
        @test isapprox(result[1, 1], 15.0; atol = 1e-10)   # t=1, z=1.5
        @test isapprox(result[1, 2], 45.0; atol = 1e-10)   # t=2, z=1.5
        @test isapprox(result[2, 1], 25.0; atol = 1e-10)   # t=1, z=2.5
        @test isapprox(result[2, 2], 55.0; atol = 1e-10)   # t=2, z=2.5
    end

    @testset "interp_along_dim: 2D (z × t) linear interpolation, 1D coord" begin
        var   = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]
        z_in  = [1.0, 2.0, 3.0]   # 1D — not same shape → interp_dim_in_is_full_array = false
        z_out = [1.5, 2.5]

        result = SSCF.interp_along_dim(
            var, 1, z_in;
            interp_dim_out = z_out,
            interp_dim_in_is_full_array = false,
        )
        @test size(result) == (2, 2)
        @test isapprox(result[1, 1], 15.0; atol = 1e-10)
        @test isapprox(result[2, 2], 55.0; atol = 1e-10)
    end

    # ------------------------------------------------------------------
    #   interp_along_dim: Missing-typed input must not crash (the key bug)
    # ------------------------------------------------------------------
    @testset "interp_along_dim: Union{Missing,Float64} input does not crash" begin
        var  = Union{Missing, Float64}[10.0 40.0; missing 50.0; 30.0 60.0]
        z_in = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out = [1.5, 2.5]

        @test_nowarn SSCF.interp_along_dim(
            var, 1, z_in;
            interp_dim_out = z_out,
            interp_dim_in_is_full_array = true,
        )
        result = SSCF.interp_along_dim(
            var, 1, z_in;
            interp_dim_out = z_out,
            interp_dim_in_is_full_array = true,
        )
        @test eltype(result) <: AbstractFloat
        @test !any(ismissing, result)
    end

    # ------------------------------------------------------------------
    #   var_to_new_coord: unweighted round-trip
    # ------------------------------------------------------------------
    @testset "var_to_new_coord: unweighted basic interpolation" begin
        var   = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]
        z_in  = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out = [1.5, 2.5]

        result = SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out)
        @test size(result) == (2, 2)
        @test isapprox(result[1, 1], 15.0; atol = 1e-10)
        @test isapprox(result[2, 2], 55.0; atol = 1e-10)
    end

    # ------------------------------------------------------------------
    #   var_to_new_coord: density-weighted gives identical result for constant weight
    # ------------------------------------------------------------------
    @testset "var_to_new_coord: uniform weight equals unweighted" begin
        var    = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]
        weight = Float64[1.0 1.0; 1.0 1.0; 1.0 1.0]   # uniform → result must equal unweighted
        z_in   = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out  = [1.5, 2.5]

        result_w  = SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out, weight)
        result_uw = SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out)
        @test isapprox(result_w, result_uw; atol = 1e-10)
    end

    # ------------------------------------------------------------------
    #   var_to_new_coord: Union{Missing,Float64} var + weight must not crash
    #   (this is the direct reproduction of the subsidence BroadcastDiskArray crash)
    # ------------------------------------------------------------------
    @testset "var_to_new_coord: Union{Missing,Float64} weighted does not crash" begin
        var    = Union{Missing, Float64}[10.0 40.0; 20.0 50.0; 30.0 60.0]
        weight = Union{Missing, Float64}[1.0 1.0; 1.0 1.0; 1.0 1.0]
        z_in   = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out  = [1.5, 2.5]

        @test_nowarn SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out, weight)
        result = SSCF.var_to_new_coord(var, z_in, 1; coord_new = z_out, weight)
        @test eltype(result) <: AbstractFloat
        @test !any(ismissing, result)
        @test isapprox(result[1, 1], 15.0; atol = 1e-10)
    end

    # ------------------------------------------------------------------
    #   var_to_new_coord: lazy-broadcast input (simulates DiskArray arithmetic)
    # ------------------------------------------------------------------
    @testset "var_to_new_coord: lazy broadcast (simulated BroadcastDiskArray) does not crash" begin
        # Simulate the production code pattern where omega has one flagged (missing) measurement:
        #   measurment at z=2, t=1 is missing  →  arithmetic keeps Union{Missing,Float64} type
        omega_nc = Union{Missing, Float64}[10.0 40.0; missing 50.0; 30.0 60.0]
        rho_nc   = Union{Missing, Float64}[1.2 1.2; 1.1 1.1; 1.0 1.0]
        # -omega / (rho * g): one missing element keeps the element type as Union{Missing,Float64}
        lazy_sub = -(omega_nc ./ (rho_nc .* 9.81))
        @test eltype(lazy_sub) >: Missing        # confirms Union{Missing,Float64}
        z_in = Float64[1.0 1.0; 2.0 2.0; 3.0 3.0]
        z_out = [1.5, 2.5]

        @test_nowarn SSCF.var_to_new_coord(lazy_sub, z_in, 1; coord_new = z_out)
        result = SSCF.var_to_new_coord(lazy_sub, z_in, 1; coord_new = z_out)
        @test eltype(result) <: AbstractFloat
        @test !any(ismissing, result)
    end

end
