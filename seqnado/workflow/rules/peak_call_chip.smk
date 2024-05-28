from typing import Literal
from seqnado.helpers import check_options
import re
import pathlib


def get_lanceotron_threshold(wildcards):
    options = config["lanceotron"]["callpeak"]
    threshold_pattern = re.compile(r"\-c\s+(\d+.?\d*)")
    threshold = threshold_pattern.search(options).group(1)
    return threshold

def format_macs_options(wildcards, options):
    query_name = f"{wildcards.sample}_{wildcards.treatment}"    
    is_paired = DESIGN.query(query_name).is_paired
    options = check_options(options)
    if not is_paired:
        options = re.sub(r"-f BAMPE", "", options)
    if not options:
        return ""
    else:
        return options


def get_control_file(wildcards, file_type: Literal["bam", "tag", "bigwig"], allow_null=False):
    control_info = DESIGN.query(sample_name=f"{wildcards.sample}_{wildcards.treatment}", full_experiment=True)["control"]
    if control_info:
        match file_type:
            case "bam":
                return f"seqnado_output/aligned/{control_info.sample_name}.bam"
            case "tag":
                return f"seqnado_output/tag_dirs/{control_info.sample_name}"
            case "bigwig":
                return f"seqnado_output/bigwigs/deeptools/unscaled/{control_info.sample_name}.bigWig"
    else:
        if allow_null:
            return []
        else:
            return "UNDEFINED"
        if allow_null:
            return []
        else:
            return "UNDEFINED"

rule macs2_with_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{treatment}.bam",
        control=lambda wc: get_control_file(wc, file_type="bam", allow_null=False),
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{treatment}.bed",
    params:
        options=lambda wc: format_macs_options(wc, config["macs"]["callpeak"]),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: f"{2 * 2 ** (attempt)}GB",
        runtime="2h",
    log:
        "seqnado_output/logs/macs/{sample}_{treatment}.log",
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -n {params.basename} {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


rule macs2_no_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{treatment}.bam",
        control=lambda wc: get_control_file(wc, file_type="bam", allow_null=True), 
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{treatment}.bed",
    params:
        options=lambda wc: format_macs_options(wc, config["macs"]["callpeak"]),
        raw=lambda wc, output: output.peaks.replace(".bed", "_peaks.xls"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: f"{2 * 2 ** (attempt)}GB",
        runtime="2h",
    log:
        "seqnado_output/logs/macs/{sample}_{treatment}.log",
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} {params.options} > {log} 2>&1 &&
        cat {params.raw} | grep -v '^#' | grep -vE '^chr\\s+start\\s+end.*' | grep -v '^$' | cut -f 1-3 > {output.peaks}
        """


rule homer_with_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}_{treatment}",
        control=lambda wc: get_control_file(wc, file_type="tag", allow_null=False),
    output:
        peaks="seqnado_output/peaks/homer/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/homer/{sample}_{treatment}.log",
    params:
        options=check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem="4GB",
        runtime="2h",
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule homer_no_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}_{treatment}",
        control=lambda wc: get_control_file(wc, file_type="tag", allow_null=True),
    output:
        peaks="seqnado_output/peaks/homer/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/homer/{sample}_{treatment}.log",
    params:
        options=check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem="4GB",
        runtime="2h",
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule lanceotron_with_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample}_{treatment}.bigWig",
        control=lambda wc: get_control_file(wc, file_type="bigwig", allow_null=False),
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/lanceotron/{sample}_{treatment}.log",
    params:
        threshold=get_lanceotron_threshold,
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    container:
        "library://asmith151/seqnado/seqnado_extra:latest"
    threads: 1
    resources:
        mem="10GB",
        runtime="6h",
    shell:
        """
        lanceotron callPeaksInput {input.treatment} -i {input.control} -f {params.outdir} --skipheader > {log} 2>&1 &&
        cat {params.basename}_L-tron.bed | awk 'BEGIN{{OFS="\\t"}} $4 >= {params.threshold} {{print $1, $2, $3}}' > {output.peaks} 
        """


rule lanceotron_no_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample}_{treatment}.bigWig",
        control=lambda wc: get_control_file(wc, file_type="bigwig", allow_null=True),
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/lanceotron/{sample}_{treatment}.bed",
    params:
        options=check_options(config["lanceotron"]["callpeak"]),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    container:
        "library://asmith151/seqnado/seqnado_extra:latest"
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt)}GB",
        runtime="6h",
    shell:
        """
        lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
        cat {params.basename}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """

rule seacr:
    input:
        treatment="seqnado_output/bedgraphs/{sample}_{treatment}.bedGraph",
    output:
        peaks="seqnado_output/peaks/seacr/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/seacr/{sample}_{treatment}.log",
    params:
        threshold=config["seacr"].get("threshold", 0.01),
        norm=config["seacr"].get("norm", "non"),
        stringency=config["seacr"].get("stringency", "stringent"),
        prefix=lambda wc, output: pathlib.Path(output.peaks).parent / pathlib.Path(output.peaks).name,
    threads: 1
    resources:
        mem="5GB",
        runtime="2h",
    shell:
        """
        SEACR_1.3.sh {input.treatment} {params.threshold} {params.norm} {params.stringency} {output.peaks} > {log} 2>&1 || touch {params.prefix}.{params.stringency}.bed
        mv {params.prefix}.{params.stringency}.bed {params.prefix}_seacr.txt
        cut -f 1-3 {params.prefix}_seacr.txt > {output.peaks}
        """
    

ruleorder: lanceotron_with_input > lanceotron_no_input
ruleorder: homer_with_input > homer_no_input
ruleorder: macs2_with_input > macs2_no_input