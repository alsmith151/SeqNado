import glob
from seqnado.helpers import get_chromosomes_from_file

# Define the list of chromosomes from your fasta index
CHROMOSOMES = get_chromosomes_from_file(config['fasta_index'])

SAMPLE_NAMES = DESIGN.sample_names

rule wasp_split_vcf:
    input:
        vcf=expand("seqnado_output/variant/bcftools/{sample}.anno.vcf.gz", sample=SAMPLE_NAMES)
    output:
        vcf=expand("seqnado_output/wasp/split_vcf/{sample}_{chromosome}.vcf.gz", sample=SAMPLE_NAMES, chromosome=CHROMOSOMES),
    wildcard_constraints:
        chromosome="chr+"
    params:
        chromosome=CHROMOSOMES, 
        prefix=expand("seqnado_output/wasp/split_vcf/{sample}", sample=SAMPLE_NAMES),
        log_prefix=expand("seqnado_output/logs/wasp/split_vcf/{sample}", sample=SAMPLE_NAMES)
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
        mem=lambda wildcards, attempt: f"{4 * 2 ** (attempt - 1)}GB"
    log:
        expand("seqnado_output/logs/wasp/split_vcf/{sample}_{chromosome}.log", sample=SAMPLE_NAMES, chromosome=CHROMOSOMES)
    shell: """
        for CHR in {params.chromosome}
        do
            bcftools view -r $CHR -Oz -o {params.prefix}_$CHR.vcf.gz {input.vcf} >> {params.log_prefix}_$CHR.log 2>&1
            bcftools index -f {params.prefix}_$CHR.vcf.gz >> {params.log_prefix}_$CHR.log 2>&1
        done
        """

rule wasp_vcf2h5:
    input:
        vcfs=expand("seqnado_output/wasp/split_vcf/{sample}_{chromosome}.vcf.gz", sample=lambda wildcards: wildcards.sample, chromosome=CHROMOSOMES)
    output:
        snp_index="seqnado_output/wasp/snp_h5/{sample}/snp_index.h5",
        snp_tab="seqnado_output/wasp/snp_h5/{sample}/snp_tab.h5",
        haplotype="seqnado_output/wasp/snp_h5/{sample}/haplotype.h5"
    wildcard_constraints:
        chromosome="chr+"
    params:
        chrom=config['fasta_index']
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
        mem=lambda wildcards, attempt: f"{4 * 2 ** (attempt - 1)}GB"
    log: "seqnado_output/logs/wasp/vcf2h5/{sample}.log"
    shell: """
        /ceph/project/milne_group/cchahrou/software/WASP/snp2h5/snp2h5 \
          --chrom {params.chrom} \
          --format vcf \
          --snp_index {output.snp_index} \
          --snp_tab {output.snp_tab} \
          --haplotype {output.haplotype} \
          {input.vcfs} 2>&1 | tee {log}
    """

rule wasp_filter_bams:
    input: rules.align_paired.output.bam
    output: "seqnado_output/wasp/mapped/{sample}.bam"
    threads: 16
    resources:
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
        mem=lambda wildcards, attempt: f"{4 * 2 ** (attempt - 1)}GB"
    log: "seqnado_output/logs/wasp/filter_bams/{sample}.log"
    shell: """
        samtools view -@ {threads} -b -q 10 {input} | \
        samtools sort -@ {threads} -o {output} >> {log} 2>&1 &&
        samtools index -@ {threads} {output} >> {log} 2>&1
        """

rule wasp_find_intersecting_snps_paired_end:
    input:
        bam=rules.wasp_filter_bams.output,
        snp_tab=rules.wasp_vcf2h5.output.snp_tab,
        snp_index=rules.wasp_vcf2h5.output.snp_index,
        haplotype=rules.wasp_vcf2h5.output.haplotype
    output:
        fastq1="seqnado_output/wasp/find_intersecting_snps/{sample}.remap.fq1.gz",
        fastq2="seqnado_output/wasp/find_intersecting_snps/{sample}.remap.fq2.gz",
        keep_bam="seqnado_output/wasp/find_intersecting_snps/{sample}.keep.bam",
        remap_bam="seqnado_output/wasp/find_intersecting_snps/{sample}.remap.bam",
    params:
        out_dir="seqnado_output/wasp/find_intersecting_snps"
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
        mem=lambda wildcards, attempt: f"{4 * 2 ** (attempt - 1)}GB",
    log: "seqnado_output/logs/wasp/find_intersecting_snps/{sample}.log"
    script: """
        ../scripts/wasp/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir {params.out_dir} \
        --snp_tab {input.snp_tab} \
        --snp_index {input.snp_index} \
        --haplotype {input.haplotype} \
        --samples {wildcards.sample} \
        {input.bam} 2> {log}
        """

use rule align_paired as wasp_remap with:
    input:
        fastq1=rules.wasp_find_intersecting_snps_paired_end.output.fastq1,
        fastq2=rules.wasp_find_intersecting_snps_paired_end.output.fastq2
    output: temp("seqnado_output/wasp/remap/{sample}.remap.bam")
    log: "seqnado_output/logs/wasp/wasp_filter_bams_remap/{sample}.log"

use rule wasp_filter_bams as wasp_filter_bams_remap with:
    input: rules.wasp_remap.output
    output: "seqnado_output/wasp/remap/{sample}.bam"
    log: "seqnado_output/logs/wasp/wasp_filter_bams_remap/{sample}.log"

rule wasp_filter_remapped_reads:
    input:
        to_remap_bam=rules.wasp_find_intersecting_snps_paired_end.output.remap_bam,
        remap_bam=rules.wasp_filter_bams_remap.output
    output:
        keep_bam="seqnado_output/wasp/filter_remapped_reads/{sample}.keep.bam"
    script: "../scripts/wasp/filter_remapped_reads.py {input.to_remap_bam} {input.remap_bam} {output.keep_bam}"

rule wasp_merge_bams:
    input:
        keep1=rules.wasp_find_intersecting_snps_paired_end.output.keep_bam,
        keep2=rules.wasp_filter_remapped_reads.output.keep_bam
    output:
        merge="seqnado_output/wasp/merge/{sample}.keep.merge.bam",
        sort="seqnado_output/wasp/merge/{sample}.keep.merge.sort.bam"
    threads: 16
    resources:
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
        mem=lambda wildcards, attempt: f"{4 * 2 ** (attempt - 1)}GB",
    log: "seqnado_output/logs/wasp/merge_bams/{sample}.log"
    shell: """
        samtools merge -@ {threads} {output.merge} {input.keep1} {input.keep2} 2> {log} &&
        samtools sort -@ {threads} -o {output.sort} {output.merge} 2>> {log} &&
        samtools index -@ {threads} {output.sort} 2>> {log}
        """

rule wasp_rmdup_pe:
    input:
        rules.wasp_merge_bams.output.sort
    output:
        rmdup="seqnado_output/wasp/rmdup/{sample}.keep.merge.rmdup.bam",
        sort="seqnado_output/wasp/rmdup/{sample}.keep.merge.rmdup.sort.bam"
    threads: 16
    resources:
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
        mem=lambda wildcards, attempt: f"{4 * 2 ** (attempt - 1)}GB",
    log: "seqnado_output/logs/wasp/rmdup_pe/{sample}.log"
    log:
    script: """
        ../scripts/wasp/rmdup_pe.py {input} {output.rmdup} 2> {log} &&
        samtools sort -@ {threads} -o {output.sort} {output.rmdup} 2>> {log} &&
        samtools index -@ {threads} {output.sort} 2>> {log}
        """

rule wasp_get_as_counts:
    input:
        bam=rules.wasp_rmdup_pe.output.sort,
        snp_index=rules.wasp_vcf2h5.output.snp_index,
        snp_tab=rules.wasp_vcf2h5.output.snp_tab,
        haplotype=rules.wasp_vcf2h5.output.haplotype
    output:
        ref_as="seqnado_output/wasp/as_counts/{sample}.ref_as_counts.h5",
        alt_as="seqnado_output/wasp/as_counts/{sample}.alt_as_counts.h5",
        other_as="seqnado_output/wasp/as_counts/{sample}.other_as_counts.h5",
        read_counts="seqnado_output/wasp/as_counts/{sample}.read_counts.h5",
        txt_counts="seqnado_output/wasp/as_counts/{sample}.as_counts.txt.gz",
    params:
        chrom=config['genome']['chromosome_sizes'],
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
        mem=lambda wildcards, attempt: f"{4 * 2 ** (attempt - 1)}GB",
    log: "seqnado_output/logs/wasp/get_as_counts/{sample}.log"
    script: """
    ../scripts/wasp/bam2h5.py \
      --chrom {params.chrom} \
      --snp_index {input.snp_index} \
      --snp_tab {input.snp_tab} \
      --ref_as_counts {output.ref_as} \
      --alt_as_counts {output.alt_as} \
      --other_as_counts {output.other_as} \
      --read_counts {output.read_counts} \
      --txt_counts {output.txt_counts} \
      {input.bam}
    """
