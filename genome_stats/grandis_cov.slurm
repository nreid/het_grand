#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

BED=~/bin/bedtools2/bin/bedtools

$BED multicov -bams $(ls /scratch/nmr15102/grandis_align/*bam | tr "\n" "\t") -bed /scratch/nmr15102/fhet_genome/5kb1kb_win.bed >/scratch/nmr15102/genome_stats/coverage/grandis_cov_5kb1kb.bed

$BED multicov -bams $(ls /scratch/nmr15102/grandis_align/*bam | tr "\n" "\t") -bed /scratch/nmr15102/fhet_genome/1kb_win.bed >/scratch/nmr15102/genome_stats/coverage/grandis_cov_1kb.bed