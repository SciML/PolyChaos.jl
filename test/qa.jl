using PolyChaos, Aqua

Aqua.find_persistent_tasks_deps(PolyChaos)
Aqua.test_ambiguities(PolyChaos, recursive = false)
Aqua.test_deps_compat(PolyChaos)
Aqua.test_piracies(PolyChaos)
Aqua.test_project_extras(PolyChaos)
Aqua.test_stale_deps(PolyChaos)
Aqua.test_unbound_args(PolyChaos, broken = true)
Aqua.test_undefined_exports(PolyChaos)
