
rule save_design:
    output:
        "seqnado_output/design.csv",
    container:
        None
    benchmark: "seqnado_output/benchmarks/save_design.benchmark",
    run:
        DESIGN.to_dataframe().to_csv("seqnado_output/design.csv", index=False)

