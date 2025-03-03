# Notes

## Notes on running the gene conversion code

Example command to run the gene conversion code:

```
configfile=config/config.yaml
threads=200
snakemake \
    --configfile $configfile \
    --cores $threads \
    --use-conda \
    -p --notemp \
    gene_conversion
```

Configuration options, special options within the config.yaml used in running the gene conversion code:

- `bed: SDs.bed`: Bed file of regions along the reference in which to try and detect gene conversion. User provided file.
- `break_paf: 10000`: This breaks the PAF on indels over 10,000 bp. You don't want significant indels to be considered as gene conversion events. Do not adjust this unless you know what you are doing.
- `min_aln_len: 1000000`: The minimum length alignment after breaking to consider when starting the gene conversion detection. Do not adjust this unless you know what you are doing.
- `window: 1000`: This is the window size used in realignment to detect gene conversion. Do not adjust this unless you know what you are doing.
- `slide: 100`: This is the slide size used in realignment to detect gene conversion. Do not adjust this unless you know what you are doing.
