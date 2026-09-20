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

    Test.@testset "ConstantBoundaryCondition is honored, not just accepted" begin
        # It is a member of `ValidBoundaryConditions`, so every signature takes it; each backend must
        # actually return the value rather than error or silently contribute nothing.
        xp = collect(0.0:1.0:4.0)
        fp = @. 1.0 + 2.0 * xp
        cbc = SSCF.Interpolation.create_bc("constant(-7.0)")
        Test.@test cbc isa SSCF.Interpolation.ConstantBoundaryCondition
        Test.@test cbc.value == -7.0

        s = SSCF.Interpolation.build_spline(
            SSCF.Interpolation.FastLinear1DInterpolation, xp, fp; bc = cbc, drop_collinear = Val(false),
        )
        Test.@test isapprox(s(2.5), 6.0; atol = 1e-12)      # in range: untouched by the bc
        Test.@test s(-1.0) == -7.0                          # below
        Test.@test s(9.0) == -7.0                           # above
        # the out-of-range part of an integral is value * width, not zero
        Test.@test isapprox(SSCF.Interpolation.safe_integrate(s, -2.0, -1.0; bc = cbc), -7.0; atol = 1e-12)
        # a single node is out of range everywhere but at the node itself
        one_node = SSCF.Interpolation.build_spline(
            SSCF.Interpolation.FastLinear1DInterpolation, [0.0], [3.0]; bc = cbc, drop_collinear = Val(false),
        )
        Test.@test one_node(0.0) == 3.0
        Test.@test one_node(5.0) == -7.0

        # the value must not widen a narrower interpolant
        s32 = SSCF.Interpolation.build_spline(
            SSCF.Interpolation.FastLinear1DInterpolation, Float32.(xp), Float32.(fp);
            bc = cbc, drop_collinear = Val(false),
        )
        Test.@test s32(2.5f0) isa Float32
        Test.@test s32(-1.0f0) isa Float32
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
