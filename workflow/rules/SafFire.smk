include: "reference_alignment.smk"


rule SafFire:
    input:
        paf=rules.trim_and_break_paf.output.paf,
    output:
        bed="results/{ref}/SafFire/{sm}.bed",
    conda:
        "../envs/env.yml"
    threads: 1
    params:
        paired_aln_len=config.get("paired_aln_len", 250_000),
        break_saffire=config.get("break_saffire", 5_000),
    shell:
        """
        rb break-paf --max-size {params.break_saffire} {input.paf} \
            | rb orient \
            | rb filter --paired-len {params.paired_aln_len} \
            | rb stats --paf  \
        > {output.bed}
        """
