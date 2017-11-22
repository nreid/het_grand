#!/bin/bash


SCRIPT1=~/het_grand/genome_stats/popcount.R
SCRIPT2=~/het_grand/genome_stats/mpile.sh

TAB=~/bin/htslib/tabix
BGZ=~/bin/htslib/bgzip
OUTDIR=/scratch/nmr15102/genome_stats/coverage_site_pop/scaffolds
OUTFILE=$1.pop.site.txt.gz

bash $SCRIPT1 $1 | Rscript $SCRIPT2 | cat | $BGZ -c >$OUTDIR/$OUTFILE