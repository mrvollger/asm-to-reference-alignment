#!/bin/bash

configfile=config/config.yaml
threads=150
snakemake --configfile $configfile --cores $threads --use-conda \
    gene_conversion \
    -p --notemp \
    "$@"

exit
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
