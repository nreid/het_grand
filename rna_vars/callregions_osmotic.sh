#!/bin/bash

scaf=$1
outdir=/scratch/nmr15102/osmotic_rnaseq/callregions

BAM=~/bin/bamtools/bin/bamtools
BED=~/bin/bedtools2/bin/bedtools
VCF=~/bin/vcflib/bin/vcfallelicprimitives


$BAM merge -list /scratch/nmr15102/osmotic_rnaseq/alignments/meta/osmotic.bams.list -region $scaf | \
samtools depth -d 200000 -q 20 -Q 30 /dev/stdin | \
awk '$3 > 1300' | awk '{OFS="\t"}{print $1,$2-1,$2,$3}' | \
$BED merge -d 20 -c 4 -o count,mean -i stdin | \
awk '$4 > 20' | \
awk '$5 < 130000' >$outdir/$scaf.bed
