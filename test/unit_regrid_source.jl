using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# Unit tests for the `AbstractRegridSource` abstraction (`AtlasInput` /
# `LESOutput`) that unifies `regrid_to_z_and_time` across data sources, and for the
# conservative regridder. Source-verb tests are pure (no NetCDF I/O); the
# conservative test uses plain arrays.
# ---------------------------------------------------------------------------
Test.@testset "Regrid source abstraction + conservative interp" begin

    Test.@testset "AbstractRegridSource: source-specific verbs dispatch per type" begin
        atlas = SSCF.AtlasInput()
        les = SSCF.LESOutput(SSCF.ObsForcing())
        Test.@test atlas isa SSCF.AbstractRegridSource
        Test.@test les isa SSCF.AbstractRegridSource

        # Atlas: pressure levels are time-varying and stored top->ground (reverse to ground->top).
        Test.@test SSCF.regrid_source_z_time_varying(atlas) == true
        Test.@test SSCF.regrid_source_reverse_z(atlas) == true
        # LES: static z column, already ground->top, starts at the initial condition.
        Test.@test SSCF.regrid_source_z_time_varying(les) == false
        Test.@test SSCF.regrid_source_reverse_z(les) == false
        Test.@test SSCF.regrid_source_initial_ind(les, nothing, 9, [0.0, 1.0, 2.0]) == 1

        # LES time is stored in days from start -> seconds from start (t - t[1]) * 86400.
        mock = Dict("time" => [10.0, 10.5, 11.0])   # days
        Test.@test SSCF.regrid_source_t_old(les, mock, Val(false)) ≈ [0.0, 0.5 * 86400, 1.0 * 86400]
        Test.@test SSCF.regrid_source_t_old(les, mock, Val(true)) == Int64[0, 86400/2, 86400]

        # data passthrough: AtlasInput returns whatever data it's given.
        d = Dict("x" => 1)
        Test.@test SSCF.regrid_source_data(atlas, 9, d) === d
    end

    Test.@testset "conservative_regridder: preserves a constant field" begin
        I = SSCF.Interpolation
        z_old = collect(0.0:100.0:1000.0)      # 11 source nodes
        f = fill(2.5, length(z_old))            # constant profile
        z_new = collect(50.0:100.0:950.0)       # 10 target nodes, within range
        out = I.conservative_regridder(z_new, z_old, f; method = I.FastLinear1DInterpolation,
            bc = I.ExtrapolateBoundaryCondition())
        Test.@test all(isapprox.(out, 2.5; atol = 1e-8))   # constant in -> constant out
    end

    Test.@testset "conservative_regridder: both modes conserve cell-integrated mass" begin
        # The regridder integrates the source spline over cells centred on the z_new nodes, with edges
        # halfway to the neighbours and extending a half-cell beyond the end nodes. Two modes:
        #   :integrate — out[i] = cell mean, so Σ out·cell_width == ∫(source spline) over the edges.
        #   :invert    — out[i] = node values whose RE-SPLINE reproduces those cell masses, so
        #                ∫(spline through (z_new, out)) over the edges == ∫(source spline) over the edges.
        # Both are tested against the same conserved total.
        # z_new is chosen INTERIOR to z_old so no cell edge extends past the data — this removes the
        # boundary-extrapolation ambiguity (safe_integrate is not additive across the data boundary),
        # making the conserved total unambiguous and forcing both modes to agree with it.
        I = SSCF.Interpolation
        z_old = collect(0.0:100.0:1000.0)
        f = 1.0 .+ 0.001 .* z_old                       # linear profile
        z_new = [250.0, 500.0, 750.0]                   # edges [125,375,625,875] ⊂ [0,1000]
        bc = I.ExtrapolateBoundaryCondition()

        # reproduce the regridder's cell edges (half-cell beyond the end nodes)
        n = length(z_new)
        edges = similar(z_new, n + 1)
        edges[1] = z_new[1] - (z_new[2] - z_new[1]) / 2
        for i in 1:(n - 1)
            edges[i + 1] = (z_new[i] + z_new[i + 1]) / 2
        end
        edges[n + 1] = z_new[n] + (z_new[n] - z_new[n - 1]) / 2
        widths = edges[2:end] .- edges[1:(end - 1)]
        @assert edges[1] ≥ z_old[1] && edges[end] ≤ z_old[end]   # no extrapolation

        spl_src = I.build_spline(I.FastLinear1DInterpolation, z_old, f; bc = bc)
        cell_means = [I.safe_integrate(spl_src, edges[i], edges[i + 1]; bc = bc) / widths[i] for i in 1:n]

        # :integrate — out[i] IS the source spline's cell mean over cell i (by construction), so
        # Σ out·width equals the summed per-cell integrals of the source spline.
        out_int = I.conservative_regridder(z_new, z_old, f; method = I.FastLinear1DInterpolation,
            bc = bc, integrate_method = :integrate)
        Test.@test isapprox(out_int, cell_means; rtol = 1e-8)
        Test.@test isapprox(sum(out_int .* widths), sum(cell_means .* widths); rtol = 1e-8)

        # :invert (default) — NOT linear-exact (A ≠ I: node values solve A·yc = cell_means); its
        # guarantee is enforce_conservation: the re-spline through (z_new, out) integrates to the same
        # total as the source spline over the edge span.
        out_inv = I.conservative_regridder(z_new, z_old, f; method = I.FastLinear1DInterpolation,
            bc = bc, integrate_method = :invert)
        respl = I.build_spline(I.FastLinear1DInterpolation, z_new, out_inv; bc = bc)
        Test.@test isapprox(I.safe_integrate(respl, edges[1], edges[end]; bc = bc),
            I.safe_integrate(spl_src, edges[1], edges[end]; bc = bc); rtol = 1e-6)
    end
end
