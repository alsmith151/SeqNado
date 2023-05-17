
use rule align_paired as align_paired_spikein with:
    params:
        index=config["genome"].get("indicies_spikein"),
        options=utils.check_options(config["bowtie2"]["options"]),
    output:
        bam=temp("seqnado_output/aligned/spikein/raw/{sample}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.log"


use rule sort_bam as sort_bam_spikein with:
    input:
        bam="seqnado_output/aligned/spikein/raw/{sample}.bam",
    output:
        bam="seqnado_output/aligned/spikein/{sample}.bam",
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.sorted.log"

use rule index_bam as index_bam_spikein with:
    input:
        bam="seqnado_output/aligned/spikein/sorted/{sample}.bam",
    output:
        bai="seqnado_output/aligned/spikein/{sample}.bam.bai",
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.sorted.index.log"



