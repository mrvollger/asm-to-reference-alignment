include: "reference_alignment.smk"


rule dip_sort:
    input:
        bam=rules.compress_sam.output.aln,
    output:
        bam=temp("temp/{ref}/bam/sorted.{sm}.bam"),
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem=8,
    shell:
        """
        samtools sort -@ {threads} -m {resources.mem} {input.bam} > {output.bam}
        """


rule dip_make_vcf:
    input:
        bam=expand(
            rules.dip_sort.output.aln,
            sm=[f"{wc.sm}_{i}" for i in [1, 2]],
            allow_missing=True,
        ),
    output:
        vcf=temp("temp/{ref}/vcf/{sm}.vcf"),
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        htsbox pileup -q0 -evcf {input.ref} {input.bam}  > {output.vcf}
        """


rule dip_phase_vcf:
    input:
        vcf=rules.dip_make_vcf.output.vcf,
    output:
        vcf="results/{ref}/vcf/{sm}.vcf.gz",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        dipcall-aux.js vcfpair -s {wildcards.sm} -a {input.vcf} | htsbox bgzip > {output.vcf}
        """


rule vcf_bed:
    input:
        vcf=rules.dip_phase_vcf.output.vcf,
    output:
        vcf="results/{ref}/vcf_bed/{sm}.bed.gz",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        #CHROM  POS0     END         ID           SVTYPE  SVLEN
        #REF    ALT     TIG_REGION  QUERY_STRAND CI      ALIGN_INDEX 
        #CLUSTER_MATCH   CALL_SOURCE     HAP     HAP_VARIANTS    GT
        bcftools norm -Ov -m-any {input.vcf} \
            | bcftools query \
                -f '%CHROM\t%POS0\t%END\t%CHROM-%POS-%TYPE-%REF-%ALT\t%TYPE\t%REF\t%ALT\t%SAMPLE\th1;h2\t%GT\n' \
                - \
            | bgzip > {output.bed}
        """
