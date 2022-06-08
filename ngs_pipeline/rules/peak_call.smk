from typing import Literal

def get_paired_ip_and_input(wc, filetype: Literal["bam", "tag", "bigwig"] = "bam"):
    
    row = (df_samples_paired.drop_duplicates(subset=["ip"])
                            .dropna()
                            .query(f"sample_name == '{wc.sample}' and antibody == '{wc.antibody}'")
                            .iloc[0]
          )

    dirs = {"bam": "aligned_and_filtered",
            "tag": "tag_dirs",
            "bigwig": "bigwigs/deeptools/"}
    
    file_ext = {"bam": ".bam",
                "tag": "/",
                "bigwig": ".bigWig"}

    try:
        return {"ip": f"{dirs[filetype]}/{row['ip']}{file_ext[filetype]}", 
                "input": f"{dirs[filetype]}/{row['input']}{file_ext[filetype]}"}
    except KeyError:
        return None


rule macs2_with_input:
    input:
        unpack(lambda wc: get_paired_ip_and_input(wc, filetype="bam")),
    output:
        peaks = "peaks/macs/{ip}.bed",
    params:
        options = config["macs"]["callpeak"],
    log:
        "logs/macs/{ip}.log",
    # conda:
    #     "../../environment_chip.yml"
    run:
        narrow = output.peaks.replace(".bed", "_peaks.narrowPeak")
        shell("""
        macs2 callpeak -t {input.ip} -c {input.input} -n peaks/macs/{wildcards.ip} -f BAM {params.options} > {log} 2>&1 &&
        cat {narrow} | cut -f 1-3 > {output.peaks}
        """)

rule macs2_no_input:
    input:
        ip = "aligned_and_filtered/{ip}.bam",
    output:
        peaks = "peaks/macs/{ip}.bed",
    params:
        options = config["macs"]["callpeak"],
    log:
        "logs/macs/{ip}.log",
    # conda:
    #     "../../environment_chip.yml"
    run:
        narrow = output.peaks.replace(".bed", "_peaks.narrowPeak")
        shell("""
        macs2 callpeak -t {input.ip} -n peaks/macs/{wildcards.ip} -f BAM {params.options} > {log} 2>&1 &&
        cat {narrow} | cut -f 1-3 > {output.peaks}
        """)


rule homer_with_input:
    input:
        unpack(lambda wc: get_paired_ip_and_input(wc, filetype="tag")),
    output:
        peaks = "peaks/homer/{sample}_{antibody}.bed",
    log:
        "logs/homer/findPeaks_{sample}_{antibody}.log",
    params:
        options = config["homer"]["findpeaks"],
    # conda:
    #     "../../environment_chip.yml"
    shell:
        """
        findPeaks {input.ip} {params.options} -o {output.peaks}.tmp  -i {input.input} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule homer_no_input:
    input:
        ip = "tag_dirs/{sample}_{antibody}/",
    output:
        peaks = "peaks/homer/{sample}_{antibody}.bed",
    log:
        "logs/homer/findPeaks_{sample}_{antibody}.log",
    params:
        options = config["homer"]["findpeaks"],
    conda:
        "../../environment_chip.yml"
    shell:
        """
        findPeaks {input.ip} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """
    

rule lanceotron_with_input:
    input:
        unpack(lambda wc: get_paired_ip_and_input(wc, filetype="bigwig")),
    output:
        peaks = "peaks/lanceotron/{ip}.bed",
    log:
        "logs/lanceotron/{ip}.log",
    params:
        options = config["lanceotron"]["callpeak"],
    shell:
        """lanceotron callPeaksInput {input.ip} -i {input.input} {params.options} -f peaks/ --format Bed --skipheader > {log} 2>&1"""

rule lanceotron_no_input:
    input:
        ip = "bigwigs/deeptools/{ip}.bigWig",
    output:
        peaks = "peaks/lanceotron/{ip}.bed",
    log:
        "logs/lanceotron/{ip}.log",
    params:
        options = config["lanceotron"]["callpeak"],
    shell:
        """lanceotron callPeaks {input.ip} {params.options} -f peaks/ --format Bed --skipheader > {log} 2>&1"""




ruleorder:
    lanceotron_with_input > lanceotron_no_input
ruleorder:
    homer_with_input > homer_no_input
ruleorder:
    macs2_with_input > macs2_no_input
