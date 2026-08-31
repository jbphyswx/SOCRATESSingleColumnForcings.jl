using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# Boundary-condition plumbing (issue #21).
#
# The regrid has two stages that ask different questions — evaluating onto `z_new` (target levels
# outside the source column) and building the returned interpolants (model times outside the forcing
# record) — so `regrid_to_z_and_time` takes `bc` as a `(; z, t)` pair. These check that the policy a
# caller asks for is the policy that actually reaches the interpolant, at the levels that can be
# exercised without artifact data.
#
# Fully qualified paths, no local aliases: test files share one `Main`, and `I` collides with
# `LinearAlgebra.I`.
# ---------------------------------------------------------------------------

Test.@testset "Boundary-condition plumbing" begin
    var = Float64[10.0 40.0; 20.0 50.0; 30.0 60.0]   # var[z, t]: 3 levels x 2 times
    z_in = [1.0, 2.0, 3.0]

    Test.@testset "built interpolants carry the requested bc" begin
        for want in (
            SSCF.Interpolation.ErrorBoundaryCondition(),
            SSCF.Interpolation.ExtrapolateBoundaryCondition(),
            SSCF.Interpolation.NearestBoundaryCondition(),
        )
            built = SSCF.interp_along_dim(
                var, 1, z_in;
                interp_dim_out = nothing,          # BUILD path: returns interpolants
                interp_dim_in_is_full_array = false,
                bc = want,
            )
            Test.@test all(itp -> itp.bc === want, built)
        end
    end

    Test.@testset "bc governs evaluation outside the node range" begin
        # 1.0 is below the first node (z_in starts at 1.0 → query 0.5 is out of range)
        nearest = SSCF.interp_along_dim(
            var, 1, z_in; interp_dim_out = [0.5], interp_dim_in_is_full_array = false,
            bc = SSCF.Interpolation.NearestBoundaryCondition(),
        )
        Test.@test isapprox(nearest[1, 1], 10.0; atol = 1e-10)   # held at the first node

        extrap = SSCF.interp_along_dim(
            var, 1, z_in; interp_dim_out = [0.5], interp_dim_in_is_full_array = false,
            bc = SSCF.Interpolation.ExtrapolateBoundaryCondition(),
        )
        Test.@test isapprox(extrap[1, 1], 5.0; atol = 1e-10)     # slope 10 per level, continued

        Test.@test_throws Exception SSCF.interp_along_dim(
            var, 1, z_in; interp_dim_out = [0.5], interp_dim_in_is_full_array = false,
            bc = SSCF.Interpolation.ErrorBoundaryCondition(),
        )
    end

    Test.@testset "var_to_new_coord forwards bc to the interpolants it builds" begin
        want = SSCF.Interpolation.NearestBoundaryCondition()
        built = SSCF.var_to_new_coord(var, z_in, 1; coord_new = nothing, bc = want)
        Test.@test all(itp -> itp.bc === want, built)
    end

    Test.@testset "conservative_mass_matrix refuses a bc it cannot honor" begin
        xc = collect(1.0:1.0:6.0)
        # The end cells extend half a spacing past the outermost centres, so the basis functions are
        # evaluated out of range by construction; an erroring bc cannot describe that.
        Test.@test_throws ErrorException SSCF.Interpolation.conservative_mass_matrix(
            xc; bc = SSCF.Interpolation.ErrorBoundaryCondition(),
        )
        for ok in (
            SSCF.Interpolation.ExtrapolateBoundaryCondition(),
            SSCF.Interpolation.NearestBoundaryCondition(),
        )
            A = SSCF.Interpolation.conservative_mass_matrix(xc; bc = ok)
            Test.@test size(A) == (length(xc), length(xc))
            Test.@test all(isfinite, A)
            Test.@test all(>=(0), A)   # clamped non-negative
        end
    end
end
