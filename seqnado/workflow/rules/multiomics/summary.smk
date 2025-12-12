
def get_assay_all_inputs():
    """Get all inputs from assay-specific 'all' rules."""
    inputs = []
    for assay in ASSAYS:
        rule_name = f"{assay}_all"
        inputs.append(getattr(rules, rule_name).input)
    return inputs


rule multiomics_summary:
    """Generate a summary report of all assays."""
    input:
        get_assay_all_inputs()
    output:
        OUTPUT_DIR + "multiomics_summary.txt"
    message:
        "Generating multiomics summary report"
    run:
        with open(output[0], 'w') as f:
            f.write("SeqNado Multiomics Project Summary\n")
            f.write("=" * 70 + "\n\n")
            f.write("Run Directory: " + str(Path.cwd()) + "\n")
            f.write(f"Total assays: {len(ASSAYS)}\n")
            f.write(f"Assays: {', '.join(ASSAYS)}\n\n")

            for assay in ASSAYS:
                f.write(f"\n{assay.upper()}\n")
                f.write("-" * 50 + "\n")
                f.write(f"  Config:   config_{assay}.yaml\n")
                f.write(f"  Metadata: metadata_{assay}.csv\n")
                f.write(f"  Output:   seqnado_output/{assay}/\n")
                f.write(f"  Report:   seqnado_output/{assay}/seqnado_report.html\n")

                # Check for completeness
                complete_file = OUTPUT_DIR + f"{assay}/seqnado_report.html"
                if os.path.exists(complete_file):
                    f.write(f"  STATUS:   COMPLETE\n")
                else:
                    f.write(f"  WARNING: SeqNado report not found! Check logs for errors.\n")
            
            f.write("\n" + "=" * 70 + "\n")
            f.write("Inputs from all assays:\n")
            for file in input:
                f.write(f"  {file}\n")
        