#!/bin/bash

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load samtools/1.3


LISTFILE=$1
OUTFILE=$(echo $1 | sed 's/.*\///' | sed 's/list/coverage.gz/' )

BGZ=~/bin/htslib/bgzip
BAM=~/bin/bamtools/bin/bamtools
BED=~/bin/bedtools2/bin/bedtools
OUTDIR=/scratch/nmr15102/genome_stats/evalsites
HICOV=/scratch/nmr15102/variants/meta/hicov.merge.sort.bed


$BAM merge -list $1 | \
$BAM filter -mapQuality ">30" -isProperPair true | \
$BED intersect -v -a stdin -b $HICOV | \
samtools depth -a /dev/stdin | \
~/bin/htslib/bgzip -c >$OUTDIR/$OUTFILE