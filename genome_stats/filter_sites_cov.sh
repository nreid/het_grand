#!/bin/bash

module load r/3.3.3

BGZ=~/bin/htslib/bgzip
BED=~/bin/bedtools2/bin/bedtools

WIN=/scratch/nmr15102/fhet_genome/1kb_win.bed
HICOV=/scratch/nmr15102/variants/meta/hicov.merge.sort.bed
COV=/scratch/nmr15102/genome_stats/coverage_site_pop/popcov.txt.gz

SCRIPT=~/het_grand/genome_stats/filter_sites_cov.R

zcat $COV | \
Rscript $SCRIPT | \
$BED intersect -v -a stdin -b $HICOV | \
$BED map -a $WIN -b stdin -c 3 -o count | \
$BGZ -c >popsites.1kb_win.bed
