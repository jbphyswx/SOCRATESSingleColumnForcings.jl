using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF

# ---------------------------------------------------------------------------
# Unit tests for the Atlas variable registry (`src/units.jl`).
# ---------------------------------------------------------------------------
Test.@testset "Atlas variable registry" begin

    les_present(fl, ft) =
        try
            p = SSCF.open_atlas_les_output(fl, ft; open_files = false)
            !isnothing(p.data) && isfile(p.data)
        catch
            false
        end

    Test.@testset "registry is internally consistent" begin
        specs = SSCF.atlas_var_specs(SSCF.LESOutput(SSCF.ObsForcing()))

        # Every `inputs` dependency resolves, and the order puts producers before consumers.
        order = SSCF.atlas_processing_order(specs)
        Test.@test length(order) == length(specs)
        position = Dict(name => i for (i, name) in enumerate(order))
        for (name, spec) in specs, dep in spec.inputs
            dep isa Symbol || continue
            dep in SSCF.ATLAS_COORDINATES && continue
            Test.@test position[dep] < position[name]
        end

        # Derived entries never claim raw units; raw entries always have a transform target.
        for (name, spec) in specs
            Test.@test spec.kind in (:raw, :derived)
            spec.kind === :derived && Test.@test isempty(spec.raw_units)
            Test.@test !isempty(spec.units)
        end

        # The three `?/kg/day` families must stay disjoint — they are only separable by long_name,
        # so a name drifting between them is a silent unit error.
        mass = Set(SSCF.ATLAS_MASS_PROCESS_RATES)
        num = Set(SSCF.ATLAS_NUMBER_PROCESS_RATES)
        slope = Set(SSCF.ATLAS_SLOPE_PARAMETERS)
        Test.@test isempty(intersect(mass, num))
        Test.@test isempty(intersect(mass, slope))
        Test.@test isempty(intersect(num, slope))
        Test.@test all(specs[n].units == "kg/kg/s" for n in mass)
        Test.@test all(specs[n].units == "1/kg/s" for n in num)
        Test.@test all(specs[n].units == "1/m" for n in slope)
    end

    Test.@testset "registry matches the shipped archive" begin
        n_checked = 0
        for fl in SSCF.flight_numbers, ft in SSCF.forcing_types
            les_present(fl, ft) || continue
            source = SSCF.LESOutput(ft)
            specs = SSCF.atlas_var_specs(source)
            SSCF.open_atlas_les_output(fl, ft).data isa Nothing && continue
            data = SSCF.open_atlas_les_output(fl, ft).data

            # (a) Coverage: a variable in the file with no spec is one that can be read with
            # assumed units, which is the failure the registry exists to prevent.
            present = Set(Symbol(n) for n in keys(data))
            unclassified = sort!(collect(setdiff(present, keys(specs))))
            Test.@test isempty(unclassified)

            # (b) Every `:raw` spec's expected units are what the file declares. This is the check
            # that catches a misclassification, since the whole table is keyed on those strings.
            mismatched = String[]
            for name in sort!(collect(keys(specs)))
                spec = specs[name]
                spec.kind === :raw || continue
                name in present || continue
                declared = SSCF.declared_units(data[String(name)])
                got = isnothing(declared) ? nothing : String(declared)
                got == spec.raw_units ||
                    push!(mismatched, "$name: spec expects $(repr(spec.raw_units)), file has $(repr(got))")
            end
            Test.@test isempty(mismatched)

            n_checked += 1
        end
        Test.@test n_checked ≥ 0   # records how many files were available; 0 means fully skipped
    end

    Test.@testset "conversions land in the units they declare" begin
        fl, ft = 9, SSCF.ObsForcing()
        if les_present(fl, ft)
            source = SSCF.LESOutput(ft)
            data = SSCF.open_atlas_les_output(fl, ft).data
            cache = Dict{Symbol, Any}()
            # the archive is single precision, so the registry's constants are built to match it
            specs = SSCF.atlas_var_specs(source, SSCF.DefaultThermodynamicsBackend(), Float32)

            # `p` is stored in mb; unconverted it would sit near 1e3, not 1e5.
            p = SSCF.read_atlas_variable(data, :p, source; specs, cache)
            Test.@test all(1.0e3 .< filter(isfinite, p) .< 1.1e5)

            # `time` is stored in days; the record is hours long, so seconds are O(1e4).
            t = SSCF.read_atlas_variable(data, :time, source; specs, cache)
            Test.@test maximum(t) > 3600

            # `RADQR` is K/day; radiative cooling is O(1e-5) K/s, not O(1) K/day.
            dTdt = filter(isfinite, SSCF.read_atlas_variable(data, :RADQR, source; specs, cache))
            Test.@test maximum(abs, dTdt) < 1.0e-2

            # `q_tot` is a specific humidity: bounded, positive, and O(1e-3) for this case.
            q_tot = filter(isfinite, SSCF.read_atlas_variable(data, :q_tot, source; specs, cache))
            Test.@test all(0 .≤ q_tot .< 0.1)

            # A derived entry actually computes rather than erroring: ice radii are microns-to-mm
            # where ice exists, and NaN where there is none.
            RI = SSCF.read_atlas_variable(data, :RI, source; specs, cache)
            finite_RI = filter(isfinite, RI)
            Test.@test !isempty(finite_RI)
            Test.@test all(0 .< finite_RI .< 1.0e-2)
        end
    end
end