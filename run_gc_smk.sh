#!/bin/bash

configfile=$2
threads=150

if [ $1 == "setup" ]; then
    snakemake --configfile $configfile --cores $threads --use-conda \
        -p --notemp \
        setup_realign "${@:3}"
fi

if [ $1 == "realign" ]; then
    snakemake --configfile $configfile --cores $threads --use-conda \
        gene_conversion \
        --config gcwindows=results/CHM13_V1.1/gene-conversion/realign/realign.bed \
        -p --notemp \
        "${@:3}"
fi
