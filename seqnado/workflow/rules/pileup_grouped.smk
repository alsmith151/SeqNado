
use rule deeptools_make_bigwigs as deeptools_make_bigwigs_consensus with:
    input:
        bam="seqnado_output/aligned/merged/{sample}.bam",
        bai="seqnado_output/aligned/merged/{sample}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/merged/{sample}.bigWig",
    threads: 8
    log:
        "seqnado_output/logs/bigwigs/{sample}.log",