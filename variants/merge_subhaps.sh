#!/bin/bash

WIN=/scratch/nmr15102/fhet_genome/100kb_win.bed
SVDIR=/scratch/nmr15102/variants/subhap/

# this line pulls out a VCF header. 
zcat /scratch/nmr15102/variants/hetgrand.vcf.gz | head -n 12000 | grep ^#

for num in $(seq 1 18873)
do
START=$(sed -n "$num"p $WIN | awk '$2=$2+1' | cut -f 2 -d ' ')
END=$(sed -n "$num"p $WIN | awk '$2=$2+1' | cut -f 3 -d ' ')
GZ=$(sed -n "$num"p $WIN | awk '$2=$2+1' | sed 's/ /_/g')
GZ=$GZ.hap.vcf.gz
zcat $SVDIR/$GZ | grep -v ^# 
done 