#!/bin/bash


SCRIPT=~/het_grand/trees_fst/fst_stream_bysite_grandis_del.R

TAB=~/bin/htslib/tabix
BGZ=~/bin/htslib/bgzip
VCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz
OUTDIR=/scratch/nmr15102/variants/fst_scaffolds_grandis_del
OUTFILE=$1.fst.gz

$TAB $VCF $1 | Rscript $SCRIPT | $BGZ -c | cat >$OUTDIR/$OUTFILE