
use rule deeptools_make_bigwigs as deeptools_make_bigwigs_consensus with:
    input:
        bam="seqnado_output/aligned/grouped/{group}.bam",
        bai="seqnado_output/aligned/grouped/{group}.bam.bai",
    output:
        bigwig="seqnado_output/bigwigs/deeptools/grouped/{group}.bigWig",
    threads: 8
    log:
        "seqnado_output/logs/bigwigs/{group}.log",