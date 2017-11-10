#!/bin/bash


SCRIPT=~/het_grand/trees_fst/fst_dxy_bysite_allpops.R

TAB=~/bin/htslib/tabix
BGZ=~/bin/htslib/bgzip
VCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz
OUTDIR=/scratch/nmr15102/variants/fst_scaffolds_allpops
OUTFILE=$1.fst.gz

$TAB $VCF $1 | Rscript $SCRIPT | $BGZ -c | cat >$OUTDIR/$OUTFILE