# Reference alignment workflow
[![DOI](https://zenodo.org/badge/414304026.svg)](https://zenodo.org/badge/latestdoi/414304026)

This repository is a snakemake workflow for aligning many genome assemblies to a reference genome using my preferred parameters, tools, and outputs. 

This workflow is also convenient for making inputs for my visualization tool [SafFire](https://mrvollger.github.io/SafFire/).

## try the test case

```
snakemake --configfile .test/config.yaml      
```

## an example run script

```
snakemake --configfile config/config.yaml 
```

## an example run script with ideograms

```
snakemake --configfile config/config.yaml ideogram 
```


### Notes on use of the pipeline in Vollger et al., 2023
Running alignment and gene conersion identification pipeline:
```
snakemake \
    --configfile config/config_asm20.yaml \
    --cores $threads \
    --use-conda \
    -p \
    gene_conversion
```
Information on where to download the input assemblies can be found on [Zenodo](https://doi.org/10.5281/zenodo.6792653).

Config files for human assemblies:
```
config/config_asm20.yaml
config/table.asm.tbl
```
Config files for the Clint PTR assembly:
```
config/clint.yaml
config/clint.asm.tbl
```
