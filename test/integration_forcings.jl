
# @warn "If you're on HPC, you should set `export OPENBLAS_NUM_THREADS=1` in your environment, or use `LinearAlgebra.BLAS.set_num_threads(1)` in your code to avoid performance issues."

# ---------------------------------------------------------------------------
# Integration tests for `get_column_forcing` and the surface functions
# (`get_surface_reference_state` / `get_surface_forcing` /
# `les_reference_profiles`). Data-guarded (skips flights/forcings whose
# input files aren't present) and NNLS-conditional (only exercises
# `enforce_positivity=true` when the NonNegLeastSquares extension is loaded).
# ---------------------------------------------------------------------------
Test.@testset "SOCRATESSingleColumnForcings integration" begin
    FT = Float64
    tp = SSCF.DefaultThermodynamicsBackend()

    # `enforce_positivity` requires the NonNegLeastSquares extension; detect it.
    _nnls_loaded =
        Base.get_extension(SSCF, :SOCRATESSingleColumnForcingsNonNegLeastSquaresExt) !== nothing
    pos_opts(cons) = cons ? (_nnls_loaded ? (false, true) : (false,)) : (false,)

    # atlas-only fields (exclude :dTdt_rad, which additionally requires an LES output file)
    atlas_vars = Tuple(v for v in SSCF.supported_forcing_variables if v != :dTdt_rad)

    files_present(fl, ft) =
        try
            p = SSCF.open_atlas_les_input(fl, ft; open_files = false)
            !isnothing(p.data) && isfile(p.data) && !isnothing(p.grid_data) && isfile(p.grid_data)
        catch
            false
        end

    new_zs = (nothing, FT.(1:100:4000))

    Test.@testset "get_column_forcing: atlas fields across the option matrix" begin
        n = 0
        for fl in SSCF.flight_numbers,
            ft in SSCF.forcing_types,
            new_z in new_zs,
            ic in (false, true),
            cons in (false, true),
            pos in pos_opts(cons)

            files_present(fl, ft) || continue
            cik = SSCF.Interpolation.get_conservative_interp_kwargs(;)
            cons && (cik = SSCF.set_property(cik, :enforce_positivity, pos))
            out = SSCF.get_column_forcing(
                fl, ft, atlas_vars, Tuple{StepRangeLen, Nothing}, Tuple{Vector, Float32};
                new_z = new_z,
                initial_condition = ic,
                thermodynamics_backend = tp,
                use_LES_output_for_z = false,
                conservative_interp = cons,
                conservative_interp_kwargs = cik,
                fail_on_missing_data = false,
            )
            Test.@test out isa NamedTuple
            Test.@test keys(out) == atlas_vars
            if !ic  # interpolant path: per-level time splines with the requested storage
                itp = out[first(atlas_vars)][1]
                Test.@test itp.xp isa AbstractRange && eltype(itp.fp) == Float32   # O(1) range time axis + Float32 value storage (specs honored)
                Test.@test isa(Test.@inferred(itp(0.0)), AbstractFloat)            # type-stable eval
            end
            n += 1
        end
        Test.@test n > 0   # at least some flight/forcing data was present and exercised
    end

    # LES output is needed for :dTdt_rad and the use_LES_output_for_z=true path.
    les_present(fl, ft) =
        try
            !isnothing(SSCF.open_atlas_les_output(fl, ft).data)
        catch
            false
        end

    Test.@testset "get_column_forcing: full field set incl. dTdt_rad — every present flight/forcing" begin
        n = 0
        for fl in SSCF.flight_numbers, ft in SSCF.forcing_types
            (files_present(fl, ft) && les_present(fl, ft)) || continue
            Test.@testset "fl=$fl $(SSCF.forcing_key(ft))" begin
                out = SSCF.get_column_forcing(fl, ft; thermodynamics_backend = tp, fail_on_missing_data = false)
                new_z = vec(Array(SSCF.default_new_z(fl)))

                # structure: one time-interpolant per requested output, at every z_new level
                Test.@test out isa NamedTuple
                Test.@test keys(out) == SSCF.supported_forcing_variables
                Test.@test isconcretetype(typeof(out))
                Test.@test all(length(getproperty(out, k)) == length(new_z) for k in keys(out))

                # physical validity of every field at the initial time
                at0(k) = [itp(0.0) for itp in getproperty(out, k)]
                for k in keys(out)
                    Test.@test all(isfinite, at0(k))
                end
                Test.@test all(250 .< at0(:H_nudge) .< 400)               # θ_liq_ice [K] over the SCM grid
                Test.@test all(180 .< at0(:T_nudge) .< 330)               # actual temperature [K]
                Test.@test all(0 .≤ at0(:qt_nudge) .< 0.05)               # total specific humidity, non-negative
                for k in (:u_nudge, :v_nudge, :ug_nudge, :vg_nudge)
                    Test.@test all(abs.(at0(k)) .< 150)                   # winds [m/s]
                end
                Test.@test all(abs.(at0(:dTdt_hadv)) .< 0.01)             # advective T tendency [K/s]
                Test.@test all(abs.(at0(:dTdt_rad)) .< 0.01)              # radiative T tendency [K/s]

                # z_old is the full-atmosphere altitude: ~0 at ground, top at tens of km, all finite.
                # ERA5 forcing is on fixed pressure levels, so the lowest 1–2 can sit below the surface
                # (p > p_sfc) → modestly negative z; the lower bound allows that but rejects unphysical values.
                zold = Array(SSCF.get_column_forcing(fl, ft; thermodynamics_backend = tp, return_old_z = true))
                Test.@test all(isfinite, zold)
                Test.@test minimum(zold) > -1e3                           # sub-surface pressure levels allowed
                Test.@test 1e4 < maximum(zold) < 6e4                      # top-of-atmosphere altitude [m]

                # use_LES_output_for_z=true derives z_old from the LES pressure profile. z_old must span all
                # forcing times, and its altitudes must be physical (0 at surface, tens of km at top, finite).
                zold_les = Array(SSCF.get_column_forcing(fl, ft; thermodynamics_backend = tp, use_LES_output_for_z = true, return_old_z = true))
                Test.@test size(zold_les, 4) == size(zold, 4)   # same time count as the hydrostatic z_old
                Test.@test all(isfinite, zold_les)
                Test.@test minimum(zold_les) > -1e3             # sub-surface pressure levels OK (matches hydrostatic)
                Test.@test 1e4 < maximum(zold_les) < 6e4
                # LES- and hydrostatic-derived altitudes agree to ~few m over the shared range (same physics)
                Test.@test isapprox(minimum(zold_les), minimum(zold); atol = 5.0)
                Test.@test !isnothing(SSCF.get_column_forcing(fl, ft; thermodynamics_backend = tp, use_LES_output_for_z = true, fail_on_missing_data = false))
            end
            n += 1
        end
        Test.@test n > 0   # at least one flight/forcing with full LES data was exercised
    end

    Test.@testset "surface reference state / conditions + LES reference profiles — every present flight/forcing" begin
        n = 0
        for fl in SSCF.flight_numbers, ft in SSCF.forcing_types
            (files_present(fl, ft) && les_present(fl, ft)) || continue
            Test.@testset "fl=$fl $(SSCF.forcing_key(ft))" begin
                # surface reference state: finite scalars; the FT arg controls the output eltype
                srs = SSCF.get_surface_reference_state(fl, ft; thermodynamics_backend = tp)
                Test.@test all(isfinite, (srs.pg, srs.Tg, srs.q_tot_g))
                Test.@test SSCF.get_surface_reference_state(fl, ft, Float32; thermodynamics_backend = tp).pg isa Float32

                # surface forcing: type-optimized interpolants — time axis stored as an O(1) range, eval finite
                # (eval alloc-freedom / inferrability of this interpolant type are enforced in allocations.jl / inferrability.jl)
                sf = SSCF.get_surface_forcing(fl, ft; thermodynamics_backend = tp)
                Test.@test keys(sf) == (:pg, :Tg, :Tsfc, :qg, :qsfc)
                Test.@test sf.pg.xp isa AbstractRange
                Test.@test all(isfinite, (sf.pg(0.0), sf.Tg(0.0), sf.Tsfc(0.0), sf.qg(0.0), sf.qsfc(0.0)))

                # LES reference profiles: structure, finiteness, positivity, eltype control, and in-place identity
                ref = SSCF.les_reference_profiles(fl; forcing_type = ft, new_zc = nothing, new_zf = nothing)
                Test.@test keys(ref) == (:p_c, :p_f, :ρ_c, :ρ_f)
                Test.@test all(all(isfinite, getproperty(ref, k)) for k in keys(ref))
                Test.@test all(all(getproperty(ref, k) .> 0) for k in keys(ref))   # p, ρ strictly positive
                Test.@test SSCF.les_reference_profiles(fl, Float32; forcing_type = ft).p_c isa Vector{Float32}
                let grid = SSCF.open_atlas_les_grid(fl).grid_data
                    nzf = [0.0; grid]; nzc = (nzf[1:(end - 1)] .+ nzf[2:end]) ./ 2
                    pc = similar(nzc); pf = similar(nzf); dc = similar(nzc); df = similar(nzf)
                    ret = SSCF.les_reference_profiles!(pc, pf, dc, df, fl; forcing_type = ft, new_zc = nzc, new_zf = nzf)
                    Test.@test ret.p_c === pc && ret.ρ_f === df   # wrote into the caller's buffers
                    Test.@test all(isfinite, pc) && all(isfinite, df)
                end
            end
            n += 1
        end
        Test.@test n > 0
    end
end
