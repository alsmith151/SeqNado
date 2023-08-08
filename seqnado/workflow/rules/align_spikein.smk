
use rule align_paired as align_paired_spikein with:
    input:
        fq1="seqnado_output/trimmed/{sample}_1.fastq.gz",
        fq2="seqnado_output/trimmed/{sample}_2.fastq.gz",
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
        bam=temp("seqnado_output/aligned/spikein/{sample, [A-Za-z0-9\-]+_[A-Za-z0-9\-]+}.bam"),
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.log"
    

use rule index_bam as index_bam_spikein with:
    input:
        bam="seqnado_output/aligned/spikein/{sample}.bam",
    output:
        bai="seqnado_output/aligned/spikein/{sample, [A-Za-z0-9\-]+_[A-Za-z0-9\-]+}.bam.bai",
    log:
        "seqnado_output/logs/aligned_spikein/{sample}.index.log"



