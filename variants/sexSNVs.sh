#!/bin/bash

source /etc/profile.d/modules.sh 

module load xz

# this script pulls out sites that distinguish males and females in heteroclitus and puts them in a file. 
# sites are those that are indentified sex-linked intervals and have a q-value < 0.001

TAB=~/bin/htslib/tabix
BGZ=~/bin/htslib/bgzip

SEXINT=/scratch/nmr15102/genome_stats/grandsex/heteroclitus_SNVsex_intervals.bed
SEXQ=/scratch/nmr15102/genome_stats/grandsex/sex_chisq_qvalues.txt.gz
VCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz
OUT=/scratch/nmr15102/variants/sexSNVs.vcf.gz

$TAB -h \
-R <($TAB -R $SEXINT $SEXQ | awk '$5 < 0.001' | awk '{OFS="\t"}{s=$2-1}{print $1,s,$2}') \
$VCF | $BGZ >$OUT

