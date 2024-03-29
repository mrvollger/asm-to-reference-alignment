window = 5000
slide = window


rule gene_conversion_target_regions:
    input:
        genome=get_fai,
        bed=config["bed"],
    output:
        bed=temp(
            expand(
                "temp/{ref}/gene-conversion/{window}-target-regions.bed",
                window=config.get("window", window),
                allow_missing=True,
            )
        ),
    threads: 1
    conda:
        "../envs/env.yml"
    params:
        slop=2 * config.get("window", window),
    shell:
        """
        bedtools sort -i {input.bed} \
            | bedtools merge -i - \
            | bedtools slop -i - -g {input.genome} -b {params.slop} \
            | bedtools merge -i - \
            > {output.bed}
        """


rule make_gene_conversion_windows:
    input:
        genome=get_fai,
        paf=ancient(rules.trim_and_break_paf.output.paf),
        bed=rules.gene_conversion_target_regions.output.bed,
    output:
        bed=temp("temp/{ref}/gene-conversion/{sm}_sliding_windows.bed"),
    log:
        "logs/{ref}/gene-conversion/{sm}_liftover.log",
    threads: 4
    conda:
        "../envs/env.yml"
    params:
        window=config.get("window", window),
        slide=config.get("slide", slide),
        min_aln_len=config.get("min_aln_len", 1e6),
        buffer=config.get("buffer", 25e3),
    shell:
        """
        awk '$4-$3>{params.min_aln_len}' {input.paf} \
            | rb --threads {threads} liftover --bed <( grep -v "^#" {input.bed}) \
            | csvtk cut  -tT -f 1,3,4 \
            | bedtools makewindows -s {params.slide} -w {params.window} -b - \
            | sed 's/$/\t{wildcards.sm}/g' \
        > {output.bed}
        """


rule make_query_windows:
    input:
        genome=get_fai,
        paf=rules.trim_and_break_paf.output.paf,
        bed=rules.make_gene_conversion_windows.output.bed,
    output:
        paf=temp("temp/{ref}/gene-conversion/{sm}_liftover.paf"),
        tbed=temp("temp/{ref}/gene-conversion/{sm}_liftover.bed"),
    log:
        "logs/{ref}/gene-conversion/{sm}_liftover.log",
    threads: 8
    conda:
        "../envs/env.yml"
    params:
        window=config.get("window", window),
    shell:
        """
        grep -w {wildcards.sm} {input.bed} \
            | cut -f 1-3 \
            | bedtools sort -i -  > {output.tbed}

        ( rb -vv --threads {threads} liftover \
            -q --bed <( grep -v "^#" {output.tbed} ) \
            --largest {input.paf} \
            | grep -v "cg:Z:{params.window}=" \
            > {output.paf} ) 2> {log}
        """


# | awk -v OFS=$'\t' '{{$2+={params.buffer}; $3-={params.buffer}}}{{print $0}}' \


rule unzip_ref:
    input:
        query=get_asm,
    output:
        fa=temp("temp/{ref}/gene-conversion/ref/{sm}_ref.fasta"),
        fai=temp("temp/{ref}/gene-conversion/ref/{sm}_ref.fasta.fai"),
    benchmark:
        "logs/{ref}/gene-conversion/alignment.{ref}_{sm}.benchmark.txt"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    shell:
        """
        seqtk seq -A -l 60 {input.query} > {output.fa}
        samtools faidx {output.fa}
        """


def get_ref(wc):
    return config.get("ref")[wc.ref]


rule window_alignment:
    input:
        #ref=rules.alignment_index.output.mmi,
        ref=get_ref,
        query=rules.unzip_ref.output.fa,
        fai=rules.unzip_ref.output.fai,
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
        minimap2 -K 1G -t {threads} \
            -cx asm20 \
            -k 15 \
            --secondary=no --eqx \
            {input.ref} \
                <( bedtools getfasta -name+ \
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
        tbl=temp("temp/{ref}/gene-conversion/{sm}_windows.tbl"),
        liftover_tbl=temp("temp/{ref}/gene-conversion/{sm}_liftover_windows.tbl"),
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        rb --threads {threads} stats --paf {input.paf} \
            > {output.tbl}
        rb --threads {threads} stats --paf {input.liftover_paf} \
            > {output.liftover_tbl}
        """


rule candidate_gene_conversion:
    input:
        window=rules.window_stats.output.tbl,
        liftover=rules.window_stats.output.liftover_tbl,
    output:
        bed=temp("temp/{ref}/gene-conversion/{sm}_candidate_windows.bed"),
    conda:
        "../envs/R.yml"
    params:
        window=config.get("window", window),
    script:
        "../scripts/combine-mappings.R"


rule group_source_windows_gene_conversion:
    input:
        bed=rules.candidate_gene_conversion.output.bed,
    output:
        bed=temp("temp/{ref}/gene-conversion/{sm}_candidate_source_windows.group.bed"),
    conda:
        "../envs/env.yml"
    params:
        find_pairs=workflow.source_path("../scripts/paired-groups.py"),
        window=config.get("window", window),
    shell:
        """
        python {params.find_pairs} \
            --fraction 0.01 --overlap 1 --source-windows \
            --distance {params.window} \
            --input {input.bed} \
        > {output.bed}
        """


rule make_query_windows_realign:
    input:
        genome=get_fai,
        paf=rules.trim_and_break_paf.output.paf,
        bed=rules.group_source_windows_gene_conversion.output.bed,
    output:
        paf=temp("temp/{ref}/gene-conversion/{sm}_liftover_2.paf"),
    log:
        "logs/{ref}/gene-conversion/{sm}_liftover_2.log",
    threads: 1
    conda:
        "../envs/env.yml"
    params:
        window=config.get("window", window),
    shell:
        """
        rb liftover -q \
            --bed <(cut -f 1-3 {input.bed} | grep -v "^#" ) \
            --largest {input.paf} \
                | grep -v "cg:Z:{params.window}=" \
        > {output.paf} 2> {log}
        """


rule window_alignment_realign:
    input:
        #ref=rules.alignment_index.output.mmi,
        ref=get_ref,
        query=rules.unzip_ref.output.fa,
        fai=rules.unzip_ref.output.fai,
        paf=rules.make_query_windows_realign.output.paf,
    output:
        aln=temp("temp/{ref}/gene-conversion/{sm}_windows_2.paf"),
    benchmark:
        "logs/{ref}/gene-conversion/alignment.{ref}_{sm}_2.benchmark.txt"
    conda:
        "../envs/env.yml"
    threads: config.get("aln_threads", 4)
    shell:
        """
        minimap2 -K 1G -t {threads} \
            -cx asm20 \
            -k 15 \
            --secondary=no --eqx \
            {input.ref} \
                <( bedtools getfasta -name+ \
                    -fi {input.query} \
                    -bed <(awk -v OFS=$'\t' '{{name=$1":"$3"-"$4}}{{print $6,$8,$9,name}}' {input.paf}) \
                ) \
            > {output.aln}
        """


rule window_stats_realign:
    input:
        paf=rules.window_alignment_realign.output.aln,
        liftover_paf=rules.make_query_windows_realign.output.paf,
    output:
        tbl=temp("temp/{ref}/gene-conversion/{sm}_windows_2.tbl"),
        liftover_tbl=temp("temp/{ref}/gene-conversion/{sm}_liftover_windows_2.tbl"),
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        rb --threads {threads} stats --paf {input.paf} \
            > {output.tbl}
        rb --threads {threads} stats --paf {input.liftover_paf} \
            > {output.liftover_tbl}
        """


rule candidate_gene_conversion_realign:
    input:
        window=rules.window_stats_realign.output.tbl,
        liftover=rules.window_stats_realign.output.liftover_tbl,
    output:
        bed=temp("temp/{ref}/gene-conversion/{sm}_candidate_windows_2.bed"),
    conda:
        "../envs/R.yml"
    params:
        window=config.get("window", window),
    script:
        "../scripts/combine-mappings.R"


rule group_gene_conversion_realign:
    input:
        bed=rules.candidate_gene_conversion_realign.output.bed,
    output:
        bed=temp("temp/{ref}/gene-conversion/{sm}_candidate_windows.group_2.bed"),
    conda:
        "../envs/env.yml"
    params:
        find_pairs=workflow.source_path("../scripts/paired-groups.py"),
    shell:
        """
        python {params.find_pairs} \
            --fraction 0.5 --overlap 1 \
            --input {input.bed} \
        > {output.bed}
        """


rule gene_conversion_windows_per_sample:
    input:
        bed=rules.group_gene_conversion_realign.output.bed,
    output:
        acceptor="results/{ref}/gene-conversion/tables/{sm}.acceptor.bed",
        bed="results/{ref}/gene-conversion/tables/{sm}.bed",
        interact="results/{ref}/gene-conversion/interactions/{sm}.bed",
    conda:
        "../envs/R.yml"
    params:
        window=config.get("window", window),
        simplify=True,
        merge=False,
    script:
        "../scripts/gene-conversion-windows.R"


rule large_table:
    input:
        beds=expand(
            rules.group_gene_conversion_realign.output.bed,
            sm=df.index,
            allow_missing=True,
        ),
    output:
        bed="results/{ref}/gene-conversion/all_candidate_windows.tbl.gz",
    conda:
        "../envs/env.yml"
    threads: 4
    shell:
        """
        (head -n 1 {input.beds[0]}; tail -q -n +2 {input.beds}) \
            | pigz -p {threads} \
            > {output.bed}
        """


rule gene_conversion_windows:
    input:
        bed=rules.large_table.output.bed,
    output:
        acceptor="results/{ref}/gene-conversion/acceptor_gene_conversion_windows.bed",
        bed="results/{ref}/gene-conversion/gene_conversion_windows.tbl",
        interact="results/{ref}/gene-conversion/gene_conversion_interactions.bed",
    conda:
        "../envs/R.yml"
    params:
        window=config.get("window", window),
        simplify=True,
        merge=False,
    script:
        "../scripts/gene-conversion-windows.R"


rule compress_large_bed:
    input:
        acceptor=rules.gene_conversion_windows.output.acceptor,
        bed=rules.gene_conversion_windows.output.bed,
        interact=rules.gene_conversion_windows.output.interact,
    output:
        acceptor=rules.gene_conversion_windows.output.acceptor + ".gz",
        bed=rules.gene_conversion_windows.output.bed + ".gz",
        interact=rules.gene_conversion_windows.output.interact + ".gz",
    conda:
        "../envs/env.yml"
    threads: 8
    shell:
        """
        bgzip -@ {threads} {input.acceptor} -c > {output.acceptor}
        bgzip -@ {threads} {input.bed} -c > {output.bed}
        bgzip -@ {threads} {input.interact} -c > {output.interact}
        """


rule make_big_bed:
    input:
        bed=rules.gene_conversion_windows.output.bed,
        interact=rules.gene_conversion_windows.output.interact,
        fai=get_fai,
    output:
        interact=(
            "results/{ref}/gene-conversion/trackHub/all_candidate_interactions.bb"
        ),
        bb="results/{ref}/gene-conversion/trackHub/all_candidate_windows.bb",
        bwa="results/{ref}/gene-conversion/trackHub/all_candidate_windows_acceptor.bw",
        bwd="results/{ref}/gene-conversion/trackHub/all_candidate_windows_donor.bw",
        bg=temp("temp/{ref}/gene-conversion/trackHub/all_candidate_windows.bg"),
        bed=temp("temp/{ref}/gene-conversion/trackHub/all_candidate_windows.bed"),
    conda:
        "../envs/env.yml"
    params:
        fmt=workflow.source_path("../scripts/bed.as"),
        interact=workflow.source_path("../scripts/interact.as"),
    threads: 1
    shell:
        """
        # make interactions
        bedtools sort -i {input.interact} \
            > {output.bed}
        bedToBigBed -as={params.interact} \
            -type=bed5+13 {output.bed} {input.fai} {output.interact}

        # all windows
        bedToBigBed -as={params.fmt} -type=bed9+10 \
            {input.bed} {input.fai} {output.bb} 

        # bed graphs
        bedtools genomecov -i <(grep -w Donor {input.bed}) \
            -g {input.fai} -bg > {output.bg}
        bedGraphToBigWig {output.bg} {input.fai} {output.bwd}

        bedtools genomecov -i <(grep -w Acceptor {input.bed}) \
            -g {input.fai} -bg > {output.bg}
        bedGraphToBigWig {output.bg} {input.fai} {output.bwa}
        """


rule make_big_beds:
    input:
        bed=rules.gene_conversion_windows_per_sample.output.bed,
        interact=rules.gene_conversion_windows_per_sample.output.interact,
        liftover_tbl=rules.window_stats_realign.output.liftover_tbl,
        fai=get_fai,
    output:
        bed=temp("temp/{ref}/gene-conversion/trackHub/gene-conversion/{sm}.bed"),
        bb="results/{ref}/gene-conversion/trackHub/gene-conversion/{sm}.bb",
        bg="results/{ref}/gene-conversion/trackHub/gene-conversion/{sm}.bg",
        interact=(
            "results/{ref}/gene-conversion/trackHub/gene-conversion/{sm}.interact.bb"
        ),
    conda:
        "../envs/env.yml"
    params:
        interact=workflow.source_path("../scripts/interact.as"),
        fmt=workflow.source_path("../scripts/bed.as"),
        slide=config.get("slide", slide),
    threads: 1
    shell:
        """
         # make interactions
         bedtools sort -i {input.interact} \
             > {output.bed}
         bedToBigBed -as={params.interact} \
             -type=bed5+13 {output.bed} {input.fai} {output.interact}

        # make others
         bedToBigBed -as={params.fmt} -type=bed9+ \
             {input.bed} {input.fai} {output.bb} 


        csvtk cut  -tT -C "$" \
            -f query_name,query_start,query_end,mismatches \
            {input.liftover_tbl} \
            | tail -n+2 \
            | awk -v OFS=$'\t' '{{print $1,$2,$2+1,$4}}' \
            | bedtools sort -i - \
            | bedtools merge -i - -c 4 -o mean \
            > {output.bed}
        bedGraphToBigWig {output.bed} {input.fai} {output.bg}
        """


rule make_trackdb:
    input:
        bigwig=expand(rules.make_big_beds.output, sm=df.index, allow_missing=True),
    output:
        track="results/{ref}/gene-conversion/trackHub/trackDb.txt",
        hub="results/{ref}/gene-conversion/trackHub/hub.txt",
        genomes="results/{ref}/gene-conversion/trackHub/genomes.txt",
    conda:
        "../envs/env.yml"
    threads: 1
    resources:
        mem=2,
        hrs=24,
    log:
        "logs/{ref}/gene-conversion/trackHub.log",
    params:
        samples=list(df.index),
    script:
        "../scripts/track_hub.py"


rule make_tbl:
    input:
        expand(
            rules.gene_conversion_windows_per_sample.output.acceptor,
            sm=df.index,
            allow_missing=True,
        ),
    output:
        tbl="results/{ref}/gene-conversion/acceptor_manifest_table.tsv",
    threads: 1
    params:
        samples=list(df.index),
    run:
        with open(output.tbl, "w+") as out:
            out.write("sample\thap\tfile\n")
            for sm, f in zip(params.samples, input):
                sm, hap = sm.split("_")
                out.write(
                    "{}\t{}\t{}\n".format(
                        sm,
                        hap,
                        os.path.abspath(f),
                    )
                )


rule merged_gene_conversion:
    input:
        bed=expand(
            rules.gene_conversion_windows_per_sample.output.acceptor,
            sm=df.index,
            allow_missing=True,
        ),
        # rules.candidate_gene_conversion_realign.output.bed,
    output:
        bed="results/{ref}/gene-conversion/merged_acceptor.bed.gz",
        tmp=temp("temp/{ref}/gene-conversion/merged_candidate_windows.group.bed.tmp"),
    conda:
        "../envs/env.yml"
    params:
        find_pairs=workflow.source_path("../scripts/paired-groups.py"),
    threads: 4
    shell:
        """
        head -n1 {input.bed[0]} > {output.tmp}
        cat {input.bed} | grep -v "^#" >> {output.tmp}

        python {params.find_pairs} \
            --fraction 0.8 --reciprocal \
            --input {output.tmp} \
            | bgzip -@ {threads} \
        > {output.bed}
        """


rule gene_conversion:
    input:
        expand(rules.make_tbl.output, ref=config.get("ref").keys()),
        expand(rules.large_table.output, ref=config.get("ref").keys()),
        expand(rules.gene_conversion_windows.output, ref=config.get("ref").keys()),
        expand(rules.compress_large_bed.output, ref=config.get("ref").keys()),
        expand(
            rules.merged_gene_conversion.output.bed,
            ref=config.get("ref").keys(),
        ),
        expand(
            rules.gene_conversion_windows_per_sample.output,
            sm=df.index,
            ref=config.get("ref").keys(),
        ),
        expand(rules.make_big_bed.output, ref=config.get("ref").keys()),
        expand(rules.make_big_beds.output, sm=df.index, ref=config.get("ref").keys()),
        expand(rules.make_trackdb.output, ref=config.get("ref").keys()),
        expand(
            rules.gene_conversion_target_regions.output,
            ref=config.get("ref").keys(),
        ),
    message:
        "Gene conversion run complete"
