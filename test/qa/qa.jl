using SciMLTesting, PolyChaos, Test

run_qa(
    PolyChaos; explicit_imports = true,
    aqua_kwargs = (; ambiguities = (; recursive = false))
)
