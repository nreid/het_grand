#!/bin/bash

# script takes one argument, a region in BED format (start 0-indexed, end 1-indexed): SCAF	START	END
# feed it a list of regions using parallel -a $LIST sh $SCRIPT
# relevant directories and files are hard-coded. 
# main purpose here is to exclude reads overlapping windows with aberrant coverage in any population.

SCAF=$(echo $1 | cut -f 1 -d ' ')
START=$(expr $(echo $1 | cut -f 2 -d ' ') + 1)
END=$(echo $1 | cut -f 3 -d ' ')

REGION=$SCAF:$START-$END
echo $REGION

OUTDIR=/scratch/nmr15102/variants/subvars
LISTFILE=/scratch/nmr15102/variants/meta/bams.list
POPSFILE=/scratch/nmr15102/variants/meta/popsfile.txt
HICOV=/scratch/nmr15102/variants/meta/hicov.merge.sort.bed

OUTFILE=$SCAF\_$START\_$END.vcf

FB=~/bin/freebayes/bin/freebayes 
BED=~/bin/bedtools2/bin/bedtools
BAM=~/bin/bamtools/bin/bamtools
REFGEN=/scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta

cd $bamdir

$BAM merge -list $LISTFILE -region $REGION | \
$BAM filter -in stdin -mapQuality ">30" -isProperPair true | \
$BED intersect -v -a stdin -b $HICOV | \
$FB -f $REFGEN --populations $POPSFILE --stdin  >$OUTDIR/$OUTFILE

echo $REGION done