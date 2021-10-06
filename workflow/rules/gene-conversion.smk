include: "reference_alignment.smk"


rule make_query_windows:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        paf="temp/gene-conversion/{ref}/{sm}.paf",
    threads: 1
    conda:
        "../envs/env.yml"
    shell:
        """
        cut -f 1,3,4 asm_small.paf \
            | bedtools makewindows -s 1000 -w 10000 -b - \
            | rustybam liftover -q --bed /dev/stdin --largest asm_small.paf \
            > {output.paf}
        """


rule window_alignment:
    input:
        ref=get_ref,
        query=get_asm,
        paf=rules.make_query_windows.output.paf,
    output:
        aln=temp("temp/gene-conversion/{ref}/{sm}.paf"),
    log:
        "logs/gene-conversion/alignment.{ref}_{sm}.log",
    benchmark:
        "logs/gene-conversion/alignment.{ref}_{sm}.benchmark.txt"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    shell:
        """
        minimap2 -K 8G -t {threads} \
            -ax asm20 \
            --secondary=no --eqx \
            {input.ref} \
                <( bedtools getfasta -fi {input.query} -bed <(cut -6,8,9 {input.paf}) ) \
            | samtools view -F 4 -b - \
            > {output.aln} 2> {log}
        """
