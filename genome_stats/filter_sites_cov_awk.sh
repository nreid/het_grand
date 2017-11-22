#!/bin/bash

module load r/3.3.3

BGZ=~/bin/htslib/bgzip
BED=~/bin/bedtools2/bin/bedtools

WIN=/scratch/nmr15102/fhet_genome/1kb_win.bed
HICOV=/scratch/nmr15102/variants/meta/hicov.merge.sort.bed
COV=/scratch/nmr15102/genome_stats/coverage_site_pop/popcov.txt.gz

SCRIPT=~/het_grand/genome_stats/filter_sites_cov.R


#min sample size: 3,3,5,5,5,3,5,5,3,5,3,3,3



zcat $COV | \
awk '	{OFS="\t"}
		{s=$2-1}
		{if(3 <  $3 && 3 < $4 && 5 < $5 && 5 < $6 && 5 < $7 && 3 < $8 && 5 < $9 && 5 < $10 && 3 < $11 && 5 < $12 && 3 < $13 && 3 < $14 && 3 < $15){print 1,$s,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}}' | \
$BED intersect -v -a stdin -b $HICOV | \
$BED map -a $WIN -b stdin -c 3 -o count | \
$BGZ -c  < popsites.1kb_win.bed.gz

