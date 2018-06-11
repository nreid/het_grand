#!/bin/bash


start=$1
stop=$(expr $start + 49)

DIR=/scratch/nmr15102/variants
OUTDIR=$DIR/rnascafvars

FB=~/bin/freebayes/bin/freebayes 
VCFI=~/bin/vcflib/bin/vcfintersect
GEN=/home/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta
BED=/scratch/nmr15102/tolerance_rnaseq/callregions_rna_union.bed
BAMLIST=/scratch/nmr15102/variants/meta/rnabams.list

$FB -k -0 --min-coverage 750 --max-coverage 300000 \
-f $GEN \
-L $BAMLIST \
-t <(sed -n $start,$stop\p $BED)  | \
$VCFI -b $BED | \
awk '$6 > 30' >$OUTDIR/bed.$start.$stop.vcf

