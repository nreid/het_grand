#!/bin/bash


BED=/home/nmr15102/bin/bedtools2/bin/bedtools
SIG=/scratch/nmr15102/genome_stats/grandsex/sex_chisq.txt.gz

zcat $SIG | \
awk '{OFS="\t"}{s=$2-1}{print $1,s,$2,$3,$4,$5,$6}' | \
awk '$4 < 0.00000001' | \
$BED merge -i stdin -d 25000 -c 4 -o count,min \
>heteroclitus_SNVsex_intervals.bed

GFF=/home/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.gff
INT=/scratch/nmr15102/genome_stats/grandsex/heteroclitus_SNVsex_intervals.bed

$BED intersect -a $GFF -b $INT | awk '$3 ~ /gene/'