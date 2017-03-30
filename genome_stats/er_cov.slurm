#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

BED=~/bin/bedtools2/bin/bedtools
BAM=~/bin/bamtools/bin//bamtools

WIN5KB=/scratch/nmr15102/fhet_genome/5kb1kb_win.bed
WIN1KB=/scratch/nmr15102/fhet_genome/1kb_win.bed

GEN=/scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.genome

LIST=/scratch/nmr15102/heteroclitus_align/meta/ER.list

OUTDIR=/scratch/nmr15102/genome_stats/coverage/

$BAM merge -list $LIST | \
$BAM filter -in - -mapQuality ">30" -isDuplicate false -isProperPair true | \
$BED bamtobed -i stdin | \
$BED map \
-a $WIN5KB \
-b stdin \
-c 5 -o mean,median,count \
-g $GEN \
>$OUTDIR/er.5kb1kb.qf.cov

$BAM merge -list $LIST | \
$BAM filter -in - -mapQuality ">30" -isDuplicate false -isProperPair true | \
$BED bamtobed -i stdin | \
$BED map \
-a $WIN1KB \
-b stdin \
-c 5 -o mean,median,count \
-g $GEN \
>$OUTDIR/er.1kb.qf.cov
