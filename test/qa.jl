using PolyChaos, Aqua, Test

Aqua.find_persistent_tasks_deps(PolyChaos)
Aqua.test_ambiguities(PolyChaos, recursive = false)
# check_extras disabled: `Pkg` is in [extras]/[targets].test without a [compat] bound.
# Tracked in https://github.com/SciML/PolyChaos.jl/issues/161
Aqua.test_deps_compat(PolyChaos, check_extras = false)
@test_broken false  # Aqua deps_compat: `Pkg` extras dep lacks [compat] — tracked in https://github.com/SciML/PolyChaos.jl/issues/161
Aqua.test_piracies(PolyChaos)
Aqua.test_project_extras(PolyChaos)
Aqua.test_stale_deps(PolyChaos)
Aqua.test_unbound_args(PolyChaos)
Aqua.test_undefined_exports(PolyChaos)
