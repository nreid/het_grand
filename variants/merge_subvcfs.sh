#!/bin/bash

WIN=/scratch/nmr15102/fhet_genome/100kb_win.bed
SVDIR=/scratch/nmr15102/variants/subvars/

zcat NW_012234324.1_1_100000.vcf.gz | grep ^#

for num in $(seq 1 18873)
do
GZ=$(sed -n "$num"p $WIN | awk '$2=$2+1' | sed 's/ /_/g')
GZ=$GZ.vcf.gz
zcat $SVDIR/$GZ | grep -v ^#
done 