#!/bin/bash



mkdir dir -p logs/drmaa

configfile=config/config.yaml
threads=200
snakemake --configfile $configfile --cores $threads --jobs $threads --use-conda \
    -p --notemp \
    dipcall gene_conversion \
    "$@"

exit
    --drmaa " -l centos=7 -l h_rt=48:00:00 -l mfree=8G -pe serial {threads} -V -cwd -S /bin/bash -w n" --drmaa-log-dir logs/drmaa \
    gene_conversion \
configfile=$2
if [ $1 == "setup" ]; then
    snakemake --configfile $configfile --cores $threads --use-conda \
        -p --notemp \
        setup_realign "${@:3}"

    cp results/CHM13_V1.1_v2.23/gene-conversion/realign/realign.bed config/premade_windows.bed
fi

if [ $1 == "realign" ]; then
    snakemake --configfile $configfile --cores $threads --use-conda \
        gene_conversion \
        --config gcwindows=config/premade_windows.bed \
        -p --notemp \
        "${@:3}"
fi
