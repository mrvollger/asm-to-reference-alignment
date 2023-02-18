# Reference alignment workflow

## an example run script

```
configfile=config/config.yaml
threads=200
snakemake --configfile $configfile --cores $threads --use-conda -p
```

Or if you want to distribute over a cluster:

```
mkdir dir -p logs/drmaa
configfile=config/config.yaml
threads=200
snakemake --configfile $configfile --jobs $threads --use-conda -p  \
    --drmaa " -l centos=7 -l h_rt=48:00:00 -l mfree=8G -pe serial {threads} -V -cwd -S /bin/bash -w n" --drmaa-log-dir logs/drmaa
```

Or if you want to make ideograms:

```
configfile=config/config.yaml
threads=200
snakemake --configfile $configfile --cores $threads --use-conda -p ideogram
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
