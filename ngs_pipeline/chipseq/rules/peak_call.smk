def get_paired_ip_and_input_bam(wc):
    row = df_samples.query(f"sample_name == '{wc.sample}' and antibody == '{wc.antibody}'").iloc[0]
    ip_file = "aligned_and_filtered/" + row["sample_name"] + "_" + row["antibody"] + ".bam"
    input_file = "aligned_and_filtered/" + row["sample_name"] + "_" + row["input"] + ".bam"
    return {"ip": ip_file, "input_control": input_file}

def get_paired_ip_and_input_tag(wc):
    row = df_samples.query(f"sample_name == '{wc.sample}' and antibody == '{wc.antibody}'").iloc[0]
    ip_file = "tag_dirs/" + row["sample_name"] + "_" + row["antibody"]
    input_file = "tag_dirs/" + row["sample_name"] + "_" + row["input"]
    return {"ip": ip_file, "input_control": input_file}


if USE_MACS:

    rule macs2_with_input:
        input:
            unpack(get_paired_ip_and_input_bam),
        output:
            peaks = "peaks/macs/{sample}_{antibody}.bed",
        params:
            options = config["macs"]["callpeak_options"] if config["macs"]["callpeak_options"] else "",
        log:
            log = "logs/macs_{sample}_{antibody}.log",
        run:
            narrow = output.peaks.replace(".bed", "_peaks.narrowPeak")
            shell("""
            macs2 callpeak -t {input.ip} -c {input.input_control} -n peaks/macs/{wildcards.sample}_{wildcards.antibody} -f BAM {params.options} 2> {log.log} &&
            cat {narrow} | cut -f 1-3 > {output.peaks}
            """)

    rule macs2_no_input:
        input:
            ip = "aligned_and_filtered/{sample}_{antibody}.bam",
        output:
            peaks = "peaks/macs/{sample}_{antibody}.bed",
        params:
            options = config["macs"]["callpeak_options"] if config["macs"]["callpeak_options"] else "",
        log:
            log = "logs/macs_{sample}_{antibody}.log",
        run:
            narrow = output.peaks.replace(".bed", "_peaks.narrowPeak")
            shell("""
            macs2 callpeak -t {input.ip} -n peaks/macs/{wildcards.sample}_{wildcards.antibody} -f BAM {params.options} 2> {log.log} &&
            cat {narrow} | cut -f 1-3 > {output.peaks}
            """)

    ruleorder: macs2_with_input > macs2_no_input >

if USE_HOMER:
    rule homer_with_input:
        input:
            ip = "tag_dirs/{sample}_{antibody}/",
            input_control = "tag_dirs/{sample}_input/",
        output:
            peaks = "peaks/homer/{sample}_{antibody}.bed",
        log:
            log = "logs/homerfindPeaks_{sample}_{antibody}.log",
        params:
            options = config["homer"]["findpeaks_options"] if config["homer"]["findpeaks_options"] else "",
        shell:
            """
            findPeaks {input.ip} {params.options} -o {output.peaks}.tmp  -i {input.input_control} 2> {log.log} &&
            pos2bed.pl {output.peaks}.tmp -o {output.peaks} &&
            rm {output.peaks}.tmp
            """


    rule homer_no_input:
        input:
            ip = "tag_dirs/{sample}_{antibody}/",
        output:
            peaks = "peaks/homer/{sample}_{antibody}.bed",
        log:
            log = "logs/homerfindPeaks_{sample}_{antibody}.log",
        params:
            options = config["homer"]["findpeaks_options"] if config["homer"]["findpeaks_options"] else "",
        shell:
            """
            findPeaks {input.ip} {params.options} -o {output.peaks}.tmp 2> {log.log} &&
            pos2bed.pl {output.peaks}.tmp -o {output.peaks} &&
            rm {output.peaks}.tmp
            """
        
    ruleorder: homer_with_input > homer_no_input 

if USE_LANCEOTRON:
    rule lanceotron_with_input:
        input:
            ip = "bigwigs/{sample}_{antibody}_deeptools.bigWig",
            input_control = "bigwigs/{sample}_input_deeptools.bigWig",
        output:
            peaks = "peaks/lanceotron/{sample}_{antibody}.bed",
        params:
            options = config["lanceotron"]["callpeak_options"] if config["lanceotron"]["callpeak_options"] else "",
        shell:
            """lanceotron callPeaksInput {input.ip} -i {input.input_control} {params.options} -f peaks/ --format Bed --skipheader """

    rule lanceotron_no_input:
        input:
            ip = "bigwigs/{sample}_{antibody}_deeptools.bigWig",
        output:
            peaks = "peaks/lanceotron/{sample}_{antibody}.bed",
        params:
            options = config["lanceotron"]["callpeak_options"] if config["lanceotron"]["callpeak_options"] else "",
        shell:
            """lanceotron callPeaks {input.ip} {params.options} -f peaks/ --format Bed --skipheader """

    ruleorder: 
        lanceotron_with_input > lanceotron_no_input