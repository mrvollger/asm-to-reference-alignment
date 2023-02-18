#!/usr/bin/env bash
set -euo pipefail
ALL_COVS=/net/eichler/vol26/projects/assembly_breaks/nobackups/asm-to-reference-alignment/human_nhp_vcf_merge/all_hap_coverages.bed                                                              
ADD_GT=/net/eichler/vol26/projects/assembly_breaks/nobackups/asm-to-reference-alignment/workflow/scripts/add-missing-hom-genotypes.py              

python $ADD_GT -e chrY -v IGC.vcf $ALL_COVS | bcftools view | bcftools sort | bgzip > IGC.fixed.gt.vcf.gz && bcftools index IGC.fixed.gt.vcf.gz 


