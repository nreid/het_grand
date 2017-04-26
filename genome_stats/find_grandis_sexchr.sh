#!/bin/bash


SCRIPT=~/het_grand/genome_stats/find_grandis_sexchr_lrat.R

TAB=~/bin/htslib/tabix
BGZ=~/bin/htslib/bgzip
VCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz
OUTDIR=/scratch/nmr15102/genome_stats/grandsex/scafsex_lrat
OUTFILE=$1.chisq.gz

$TAB $VCF $1 | Rscript $SCRIPT | $BGZ -c | cat >$OUTDIR/$OUTFILE