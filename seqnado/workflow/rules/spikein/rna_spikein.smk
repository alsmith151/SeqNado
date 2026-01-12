from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested

use rule align_paired as align_paired_spikein_rna with:
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.sortedByCoord.out.bam"),
        bam2=temp(
            OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.toTranscriptome.out.bam"
        ),
        log_out=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Log.final.out"),
    params:
        index=str(CONFIG.genome.index.prefix),
        options=str(CONFIG.third_party_tools.star.align.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/spikein/{sample}_",
    log: OUTPUT_DIR + "/logs/align_spikein/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align_spikein/{sample}.tsv",

use rule align_single as align_single_spikein_rna with:
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.sortedByCoord.out.bam"),
        bam2=temp(
            OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.toTranscriptome.out.bam"
        ),
        log_out=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Log.final.out"),
    params:
        index=str(CONFIG.genome.index.prefix),
        options=str(CONFIG.third_party_tools.star.align.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/spikein/{sample}_",
    log: OUTPUT_DIR + "/logs/align_spikein/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align_spikein/{sample}.tsv",

use rule rename_aligned as rename_aligned_spikein with:
    input:
        bam=OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/raw/{sample}.bam"),

localrules:
    rename_aligned_spikein,
