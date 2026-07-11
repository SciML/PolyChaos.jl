using SciMLTesting, PolyChaos, Test

run_qa(
    PolyChaos; explicit_imports = true,
    api_docs_kwargs = (; rendered = true),
    aqua_kwargs = (; ambiguities = (; recursive = false))
)
