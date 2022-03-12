include: "reference_alignment.smk"


rule pav_bam:
    input:
        paf=rules.trim_and_break_paf.output.paf,
    output:
        bam="results/{ref}/pav_input_bam/{sm}.bam",
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem=8,
    params:
        min_aln_len=config.get("min_aln_len", 1e6),
    shell:
        """
        awk '$4-$3>{params.min_aln_len}' {input.paf} \
            | rb paf-to-sam \
            | samtools sort -@ {threads} -m {resources.mem}G \
        > {output.bam}
        """


rule dip_sort:
    input:
        paf=rules.trim_and_break_paf.output.paf,
        query=get_asm,
    output:
        bam=temp("temp/{ref}/bam/sorted.{sm}.bam"),
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem=8,
    params:
        min_aln_len=config.get("min_aln_len", 1e6),
    shell:
        """
        awk '$4-$3>{params.min_aln_len}' {input.paf} \
            | rb paf-to-sam -f {input.query} \
            | samtools sort -@ {threads} -m {resources.mem}G \
        > {output.bam}
        """


# used as a test only
rule dip_sort_bam:
    input:
        aln=rules.compress_sam.output.aln,
    output:
        bam=temp("temp/{ref}/bam/bam.sorted.{sm}.bam"),
    conda:
        "../envs/env.yml"
    threads: 4
    resources:
        mem=8,
    shell:
        """
        samtools sort -@ {threads} -m {resources.mem}G {input.aln} > {output.bam}
        """


rule dip_make_vcf:
    input:
        bam=lambda wc: expand(
            rules.dip_sort.output.bam,
            sm=[f"{wc.sm}_{i}" for i in [1, 2]],
            allow_missing=True,
        ),
        ref=get_ref,
    output:
        vcf=temp("temp/{ref}/vcf/{sm}.vcf"),
    conda:
        "../envs/dipcall.yml"
    threads: 1
    shell:
        """
        htsbox pileup -q0 -evcf {input.ref} {input.bam}  > {output.vcf}
        """


rule dip_phase_vcf:
    input:
        vcf=rules.dip_make_vcf.output.vcf,
        ref=get_ref,
    output:
        vcf="results/{ref}/vcf/{sm}.vcf.gz",
    conda:
        "../envs/dipcall.yml"
    log:
        "logs/{ref}/dipcall_phased_vcf/{sm}.log",
    resources:
        mem=8,
    threads: 1
    shell:
        """
        ( dipcall-aux.js vcfpair -s {wildcards.sm} -a {input.vcf} \
            | bcftools norm -Ov -m-any \
            | bcftools norm -Ov -d exact \
            | bcftools norm -Ov -m-any --fasta-ref {input.ref} --check-ref w \
            | bcftools sort -m {resources.mem}G \
            | htsbox bgzip > {output.vcf} ) 2> {log}
        """


rule dip_phase_index:
    input:
        vcf=rules.dip_phase_vcf.output.vcf,
    output:
        idx=f"{rules.dip_phase_vcf.output.vcf}.csi",
    conda:
        "../envs/dipcall.yml"
    shell:
        """
        bcftools index {input.vcf}
        """


rule merged_vcf:
    input:
        vcf=expand(
            rules.dip_phase_vcf.output.vcf,
            sm=df["sample"].str.strip().unique(),
            ref=config.get("ref").keys(),
        ),
        idx=expand(
            rules.dip_phase_index.output,
            sm=df["sample"].str.strip().unique(),
            ref=config.get("ref").keys(),
        ),
    output:
        vcf="results/{ref}/merged.vcf.gz",
    conda:
        "../envs/dipcall.yml"
    threads: 8
    params:
        n_samples=len(df["sample"].str.strip().unique()),
    shell:
        """
        if [ {params.n_samples} == 1 ]; then
            bcftools norm --threads {threads} -Ov -m-any {input.vcf} \
                | bgzip -@ {threads} \
                > {output.vcf}
        else 
            bcftools merge \
                --threads {threads} {input.vcf} \
                | bcftools norm --threads {threads} -Ov -m-any \
                | bgzip -@ {threads} \
                > {output.vcf}
        fi
        """


rule callable_regions:
    input:
        paf=rules.trim_and_break_paf.output.paf,
    output:
        bed="results/{ref}/callable/{sm}_callable_regions.bed.gz",
    conda:
        "../envs/env.yml"
    threads: 1
    params:
        min_aln_len=config.get("min_aln_len", 1e6),
    shell:
        """
        awk '$4-$3>{params.min_aln_len}' {input.paf} \
            | csvtk cut  -tT -f 6,8,9,1,3,4 \
            | bgzip > {output.bed}
        """


rule add_missing_genotypes:
    input:
        vcf=rules.merged_vcf.output.vcf,
        bed=expand(
            rules.callable_regions.output,
            sm=df.index,
            ref=config.get("ref").keys(),
        ),
    output:
        bed=temp("results/{ref}/callable/callable_regions.bed"),
        vcf="results/{ref}/merged.hom.fixed.vcf.gz",
    conda:
        "../envs/env.yml"
    threads: 1
    params:
        add_gt=workflow.source_path("../scripts/add-missing-hom-genotypes.py"),
    shell:
        """
        ls {input.bed} \
            | parallel -n 1 \
            $'zcat {} | cut -f 1-3 | sed \'s/$/\\t{{/.}}/g\' '
        > {output.bed}
        head {output.bed}

        python {params.add_gt} -v {input.vcf} {output.bed} \
            | bgzip -@ {threads} \
            > {output.vcf}
        """


rule vcf_bed:
    input:
        vcf=rules.dip_phase_vcf.output.vcf,
        ref=get_ref,
    output:
        bed="results/{ref}/vcf_bed/{sm}.bed.gz",
    conda:
        "../envs/env.yml"
    threads: 1
    params:
        header="#CHROM\tPOS\tEND\tID\tTYPE\tREF\tALT\tSAMPLE\tHAP\tGT",
    shell:
        """
        ( echo '{params.header}'; \
            bcftools query \
                    -f '%CHROM\t%POS0\t%END\t%CHROM-%POS-%TYPE-%REF-%ALT\t%TYPE\t%REF\t%ALT\t{wildcards.sm}\th1;h2\t[ %GT]\n' \
                    {input.vcf} \
        ) \
            | bgzip > {output.bed}
        """
