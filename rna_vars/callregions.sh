#!/bin/bash

scaf=$1
outdir=/scratch/nmr15102/tolerance_rnaseq/callregions

BAM=~/bin/bamtools/bin/bamtools
BED=~/bin/bedtools2/bin/bedtools
VCF=~/bin/vcflib/bin/vcfallelicprimitives


$BAM merge -list /scratch/nmr15102/tolerance_rnaseq/alignments/meta/bams.list -region $scaf | \
samtools depth -d 200000 -q 20 -Q 30 /dev/stdin | \
awk '$3 > 750' | awk '{OFS="\t"}{print $1,$2-1,$2,$3}' | \
$BED merge -d 20 -c 4 -o count,mean -i stdin | \
awk '$4 > 20' | \
awk '$5 < 75000' >$outdir/$scaf.bed
