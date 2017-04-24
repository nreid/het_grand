#!/bin/bash


SCRIPT=~/het_grand/trees_fst/fst_dxy_bysite.R

TAB=~/bin/htslib/tabix
BGZ=~/bin/htslib/bgzip
VCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz
OUTDIR=/scratch/nmr15102/variants
OUTFILE=$1.fst.gz

$TAB $VCF $1 | Rscript $SCRIPT | $BGZ -c >$OUTDIR/$OUTFILE