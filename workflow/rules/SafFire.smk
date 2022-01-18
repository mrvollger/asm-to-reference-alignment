include: "reference_alignment.smk"


rule SafFire:
    input:
        paf=rules.sam_to_paf.output.paf,
    output:
        bed="results/{ref}/SafFire/{sm}.bed",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        rb break-paf --max-size 5000 {input.paf} \
            | rb orient \
            | rb filter --paired-len 250000 \
            | rb stats --paf  \
        > {output.bed}
        """
