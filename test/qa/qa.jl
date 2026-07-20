using SciMLTesting, PolyChaos, Test

run_qa(
    PolyChaos;
    aqua_kwargs = (; ambiguities = (; recursive = false))
)
