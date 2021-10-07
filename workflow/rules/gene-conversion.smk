include: "reference_alignment.smk"


rule make_query_windows:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        paf=temp("temp/{ref}/gene-conversion/{sm}_liftover.paf"),
    threads: 1
    conda:
        "../envs/env.yml"
    params:
        window=config.get("window", 10000),
        slide=config.get("slide", 1000),
    shell:
        """
        cut -f 1,3,4 {input.paf} \
            | bedtools makewindows -s {params.slide} -w {params.window} -b - \
            | rb liftover -q --bed /dev/stdin --largest {input.paf} \
            > {output.paf}
        """


rule unzip_ref:
    input:
        query=get_asm,
    output:
        temp("temp/{ref}/gene-conversion/ref/{sm}_ref.fasta"),
    benchmark:
        "logs/{ref}/gene-conversion/alignment.{ref}_{sm}.benchmark.txt"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    shell:
        """
        seqtk seq -A -l 60 {input.query} > {output}
        samtools index {output}
        """


rule window_alignment:
    input:
        ref=get_ref,
        query=rules.unzip_ref.output,
        paf=rules.make_query_windows.output.paf,
    output:
        aln=temp("temp/{ref}/gene-conversion/{sm}_windows.paf"),
    benchmark:
        "logs/{ref}/gene-conversion/alignment.{ref}_{sm}.benchmark.txt"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    shell:
        """
        minimap2 -K 8G -t {threads} \
            -cx asm20 \
            --secondary=no --eqx \
            {input.ref} \
                <( bedtools getfasta -name \
                    -fi {input.query} \
                    -bed <(awk -v OFS=$'\t' '{{name=$1":"$3"-"$4}}{{print $6,$8,$9,name}}' {input.paf}) \
                ) \
            > {output.aln}
        """


rule window_stats:
    input:
        paf=rules.window_alignment.output.aln,
        liftover_paf=rules.make_query_windows.output.paf,
    output:
        tbl="results/{ref}/gene-conversion/{sm}_windows.tbl.gz",
        liftover_tbl="results/{ref}/gene-conversion/{sm}_liftover_windows.tbl.gz",
    conda:
        "../envs/env.yml"
    threads: 4
    shell:
        """
        rb stats --paf {input.paf} \
            | pigz -p 4  > {output.tbl}
        rb stats --paf {input.liftover_paf} \
            | pigz -p {threads} > {output.liftover_tbl}
        """


rule candidate_gene_conversion:
    input:
        window=rules.window_stats.output.tbl,
        liftover=rules.window_stats.output.liftover_tbl,
    output:
        tbl="results/{ref}/gene-conversion/{sm}_candidate_windows.tbl",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/combine-mappings.R"


rule gene_conversion:
    input:
        expand(
            rules.candidate_gene_conversion.output,
            sm=df.index,
            ref=config.get("ref").keys(),
        ),
    message:
        "Gene conversion run complete"
