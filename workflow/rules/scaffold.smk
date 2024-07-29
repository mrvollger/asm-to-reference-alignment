rule run_scaffold:
    input:
        paf=rules.sam_to_paf.output.paf,
        ref=get_ref,
        asm=get_asm,
    output:
        fa="results/{ref}/scaffold/{sm}.scaffold.fa",
    shadow:
        "minimal"
    conda:
        "../envs/ragtag.yaml"
    params:
        g=config.get("ragtag_g", 100),
        m=config.get("ragtag_m", 10_000_000),
        aln_len=config.get("ragtag_aln_len", 1_000_000),
    shell:
        """
        # make a directory for the ragtag output
        mkdir -p ragtag

        # add the paf we want to use for scaffolding
        rb filter -p {params.aln_len} {input.paf} > ragtag/ragtag.scaffold.asm.paf 

        # run ragtag
        ragtag.py scaffold \
            -u -r -g {params.g} -m {params.m} \
            -o ragtag \
            {input.ref} {input.asm} 
        
        # add in the sample tag
        sed -i 's/_RagTag/_{wildcards.sm}/' ragtag/ragtag.scaffold.fasta

        # move the scaffolded fasta to the output 
        mv ragtag/ragtag.scaffold.fasta {output.fa}
        """


rule index_scaffold:
    input:
        fa=rules.run_scaffold.output.fa,
    output:
        fai=f"{rules.run_scaffold.output.fa}.fai",
    conda:
        "../envs/env.yml"
    shell:
        """
        samtools faidx {input.fa}
        """


rule scaffold_gaps:
    input:
        fa=rules.run_scaffold.output.fa,
    output:
        bed="results/{ref}/scaffold/{sm}_gaps.bed",
    conda:
        "../envs/env.yml"
    shell:
        """
        printf "#ct\tst\ten\tname\tscore\tstrand\ttst\tten\tcolor\tb_ct\tb_sz\tb_st\n" > {output.bed}
        seqtk gap {input.fa} >> {output.bed}
        """


rule scaffold:
    input:
        fai=expand(
            rules.index_scaffold.output.fai,
            ref=config.get("ref").keys(),
            sm=df.index,
        ),
        fa=expand(
            rules.run_scaffold.output.fa,
            ref=config.get("ref").keys(),
            sm=df.index,
        ),
        gaps=expand(
            rules.scaffold_gaps.output.bed,
            ref=config.get("ref").keys(),
            sm=df.index,
        ),
    shell:
        """
        for fai in {input.fai}; do
            echo $fai
            grep -v '^h' $fai | awk '{{sum+=$2}} END {{print sum/1e6}}'
            grep '^h' $fai | awk '{{sum+=$2}} END {{print sum/1e6}}'
        done

        echo ""
        
        for gap in {input.gaps}; do
            rb bl -r $gap
        done
        """
