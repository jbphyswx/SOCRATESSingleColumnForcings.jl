using Test
import ClimaParams as CP
import SOCRATESSingleColumnForcings as SSCF
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

# ---------------------------------------------------------------------------
# Running policy
#
#   Default (no args):        run ALL tests — unit and integration.
#   Opt-out from integration: pass "no_integration" as a test arg, OR set
#                             SSCF_SKIP_INTEGRATION_TESTS=true in the env.
#   Run only unit tests:      pass "unit" as a test arg.
#   Run only integration:     pass "integration" as a test arg.
#   Run specific unit group:  pass "unit_ncdatasets", "unit_shapes", or
#                             "unit_interp" as a test arg.
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
    include("unit_ncdatasets_compat.jl")
end
if _run_unit_group("unit_shapes")
    include("unit_shape_contracts.jl")
end
if _run_unit_group("unit_interp")
    include("unit_interp_methods.jl")
end

if !_skip_integration()
    include("integration_process_case.jl")
end
