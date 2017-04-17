#!/bin/bash

#module load r/3.2.3
# NW_012224401.1

TAB=~/bin/htslib/tabix

HAPVCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz

DIR=/scratch/nmr15102/sweepfinder
OUT=fst.null.out

$TAB $HAPVCF $1 | Rscript ~/het_grand/trees_fst/fst_stream.R >$DIR/$OUT