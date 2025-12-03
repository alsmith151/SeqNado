match ASSAY:
    case Assay.RNA:
        include: "rna_spikein.smk"
    case Assay.ATAC | Assay.CAT | Assay.CHIP:
        include: "dna_spikein.smk"
    case _:
        raise ValueError(f"Spike-in not supported for assay {ASSAY}")

include: "spikein.smk"
