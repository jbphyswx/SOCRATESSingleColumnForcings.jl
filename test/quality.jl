using Test: Test
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
import Aqua
import ExplicitImports

# ---------------------------------------------------------------------------
# Package-quality gates: Aqua (ambiguities, unbound args, stale/undeclared deps, compat bounds,
# piracy, project extras) and ExplicitImports (no implicit imports, no stale explicit imports).
# These guard the packaging invariants that behavior tests cannot see.
# ---------------------------------------------------------------------------
Test.@testset "Aqua" begin
    Aqua.test_ambiguities(SSCF)
    Aqua.test_unbound_args(SSCF)
    Aqua.test_undefined_exports(SSCF)
    Aqua.test_stale_deps(SSCF)
    Aqua.test_deps_compat(SSCF)
    Aqua.test_piracies(SSCF)
    Aqua.test_project_extras(SSCF)
    Aqua.test_persistent_tasks(SSCF)
end

Test.@testset "ExplicitImports" begin
    Test.@test ExplicitImports.check_no_implicit_imports(SSCF) === nothing
    Test.@test ExplicitImports.check_no_stale_explicit_imports(SSCF) === nothing
end
