using Test: Test
using ClimaParams: ClimaParams as CP
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Thermodynamics: Thermodynamics as TD
using NonNegLeastSquares: NonNegLeastSquares  # activates the NNLS extension so enforce_positivity is exercised

# ---------------------------------------------------------------------------
# Running policy
#
#   Default (no args):        run ALL tests — unit and integration.
#   Opt-out from integration: pass "no_integration" as a test arg, OR set
#                             SSCF_SKIP_INTEGRATION_TESTS=true in the env.
#   Run only unit tests:      pass "unit" as a test arg.
#   Run only integration:     pass "integration" as a test arg.
#   Run specific unit group:  pass "unit_ncdatasets", "unit_shapes",
#                             "unit_interp", "unit_thermo", or "unit_regrid".
# ---------------------------------------------------------------------------
const _test_args = Set(ARGS)

_skip_integration() =
    ("no_integration" in _test_args) ||
    get(ENV, "SSCF_SKIP_INTEGRATION_TESTS", "false") == "true"

# True when the user asked for *only* integration (no unit flag)
_only_integration() =
    ("integration" in _test_args) &&
    !("unit" in _test_args) &&
    !any(a -> startswith(a, "unit_"), _test_args)

_run_unit_group(name::AbstractString) =
    !_only_integration() &&
    (isempty(_test_args) || "unit" in _test_args || name in _test_args)

if _run_unit_group("unit_ncdatasets")
    include("unit_ncdatasets.jl")
end
if _run_unit_group("unit_shapes")
    include("unit_shape_contracts.jl")
end
if _run_unit_group("unit_interp")
    include("unit_interp_methods.jl")
end
if _run_unit_group("unit_thermo")
    include("unit_thermodynamics.jl")
end
if _run_unit_group("unit_thermo_consistency")
    include("unit_thermodynamics_consistency.jl")
end
if _run_unit_group("unit_regrid")
    include("unit_regrid_source.jl")
end
if _run_unit_group("quality")
    include("quality.jl")
end
if _run_unit_group("extensions")
    include("extensions.jl")
end
if _run_unit_group("inferrability")
    include("inferrability.jl")
end
if _run_unit_group("allocations")
    include("allocations.jl")
end

if !_skip_integration()
    include("integration_forcings.jl")
end
