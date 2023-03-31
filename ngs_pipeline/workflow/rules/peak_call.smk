from typing import Literal


def get_paired_treatment_and_input(
    wc, filetype: Literal["bam", "tag", "bigwig"] = "bam"
):

    paths = dict()

    dirs = {
        "bam": "aligned_and_filtered",
        "tag": "tag_dirs",
        "bigwig": "bigwigs/deeptools/",
    }

    file_ext = {"bam": ".bam", "tag": "/", "bigwig": ".bigWig"}

    paths["treatment"] = os.path.join(
        dirs[filetype], f"{wc.treatment}{file_ext[filetype]}"
    )

    if (ASSAY == "ChIP") and (
        config["lanceotron"].get("use_input", True)
        or config["homer"].get("use_input", True)
    ):

        if SAMPLE_NAMES_PAIRED[wc.treatment]:
            paths["control"] = os.path.join(
                dirs[filetype],
                f"{SAMPLE_NAMES_PAIRED[wc.treatment]}{file_ext[filetype]}",
            )

    else:
        paths["control"] = ""

    return paths


rule macs2_with_input:
    input:
        unpack(lambda wc: get_paired_treatment_and_input(wc, filetype="bam")),
    output:
        peaks="peaks/macs/{treatment}.bed",
    params:
        options=config["macs"]["callpeak"],
        narrow=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
    log:
        "logs/macs/{treatment}.log",
    shell:
        # # narrow = output.peaks.replace(".bed", "_peaks.narrowPeak")
        # params.narrow
        # cmd = f"""
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -n peaks/macs/{wildcards.treatment} -f BAM {params.options} > {log} 2>&1 &&
        cat {params.narrow} | cut -f 1-3 > {output.peaks}
        """
        # """
        # if workflow.use_singularity:
        #         cmd = utils.get_singularity_command(command=cmd,
        #                                             workflow=workflow,)
        # shell(cmd)


rule macs2_no_input:
    input:
        treatment="aligned_and_filtered/{treatment}.bam",
    output:
        peaks="peaks/macs/{treatment}.bed",
    params:
        options=config["macs"]["callpeak"],
    log:
        "logs/macs/{treatment}.log",
    run:
        narrow = output.peaks.replace(".bed", "_peaks.narrowPeak")
        cmd = f"""
                                                macs2 callpeak -t {input.treatment} -n peaks/macs/{wildcards.treatment} -f BAM {params.options} > {log} 2>&1 &&
                                                cat {narrow} | cut -f 1-3 > {output.peaks}
                                                """

        if workflow.use_singularity:
            cmd = utils.get_singularity_command(
                command=cmd,
                workflow=workflow,
            )
        shell(cmd)


rule homer_with_input:
    input:
        unpack(lambda wc: get_paired_treatment_and_input(wc, filetype="tag")),
    output:
        peaks="peaks/homer/{treatment}.bed",
    log:
        "logs/homer/findPeaks/{treatment}.log",
    params:
        options=config["homer"]["findpeaks"],
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule homer_no_input:
    input:
        treatment="tag_dirs/{sample}_{antibody}/",
    output:
        peaks="peaks/homer/{sample}_{antibody}.bed",
    log:
        "logs/homer/findPeaks_{sample}_{antibody}.log",
    params:
        options=str(config["homer"]["findpeaks"]).replace("None", ""),
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule lanceotron_with_input:
    input:
        unpack(lambda wc: get_paired_treatment_and_input(wc, filetype="bigwig")),
    output:
        peaks="peaks/lanceotron/{treatment}.bed",
    log:
        "logs/lanceotron/{treatment}.log",
    params:
        options=config["lanceotron"]["callpeak"],
    threads: 1
    resources:
        mem_mb=1024 * 10,
    shell:
        """lanceotron callPeaksInput {input.treatment} -i {input.control} -f peaks/ --skipheader > {log} 2>&1 &&
           cat peaks/{wildcards.treatment}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """


rule lanceotron_no_input:
    input:
        treatment="bigwigs/deeptools/{treatment}.bigWig",
    output:
        peaks="peaks/lanceotron/{treatment}.bed",
    log:
        "logs/lanceotron/{treatment}.log",
    params:
        options=config["lanceotron"]["callpeak"],
    resources:
        mem_mb=1024 * 10,
    shell:
        """lanceotron callPeaks {input.treatment} {params.options} -f peaks/ --format Bed --skipheader > {log} 2>&1 &&
           mv peaks/{wildcards.treatment}_L-tron.bed {output.peaks}
        """


ruleorder: lanceotron_with_input > lanceotron_no_input
ruleorder: homer_with_input > homer_no_input
ruleorder: macs2_with_input > macs2_no_input
