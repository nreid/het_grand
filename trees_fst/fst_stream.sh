#!/bin/bash

#module load r/3.2.3
# NW_012224401.1

TAB=~/bin/htslib/tabix

HAPVCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz

$TAB $HAPVCF $1 | Rscript ~/het_grand/trees_fst/fst_stream.R 