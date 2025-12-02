from seqnado.helpers import define_memory_requested, define_time_requested

use rule align_paired as align_paired_spikein with:
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.sortedByCoord.out.bam"),
        bam2=temp(
            OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.toTranscriptome.out.bam"
        ),
        log_out=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Log.final.out"),

use rule align_single as align_single_spikein with:
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.sortedByCoord.out.bam"),
        bam2=temp(
            OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.toTranscriptome.out.bam"
        ),
        log_out=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Log.final.out"),

use rule rename_aligned as rename_aligned_spikein with:
    input:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_Aligned.sortedByCoord.out.bam"),
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/raw/{sample}.bam"),

localrules:
    rename_aligned_spikein,
