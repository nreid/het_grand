#!/bin/bash


# This script executes an R script that haploidizes our VCF file genotypes. 
 	# zcat hetgrand.vcf.gz | grep -v ^# | Rscript test.R >test.out

SCAF=$(echo $1 | cut -f 1 -d ' ')
START=$(expr $(echo $1 | cut -f 2 -d ' ') + 1)
END=$(echo $1 | cut -f 3 -d ' ')

REGION=$SCAF:$START-$END
echo $REGION

HAP=~/het_grand/variants/haploidize_stream.R
OUTDIR=/scratch/nmr15102/variants/subhap

VCF=/scratch/nmr15102/variants/hetgrand.vcf.gz 

OUTFILE=$SCAF\_$START\_$END.hap.vcf
echo $OUTFILE

TAB=~/bin/htslib/bgzip
BGZ=~/bin/htslib/bgzip


$TAB $VCF $REGION | Rscript $HAP >$OUTDIR/$OUTFILE

$BGZ $OUTDIR/$OUTFILE


