
rule save_design:
    output:
        "seqnado_output/design.csv",
    container:
        None
    benchmark: "seqnado_output/benchmark/design.benchmark" if config.get("benchmark", False) else None
    run:
        DESIGN.to_dataframe().to_csv("seqnado_output/design.csv", index=False)

