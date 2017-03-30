#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

BED=~/bin/bedtools2/bin/bedtools

$BED nuc -fi /scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta -bed /scratch/nmr15102/fhet_genome/5kb1kb_win.bed >/scratch/nmr15102/genome_stats/basecomp/ACGT_5kb1kb.bed

$BED nuc -fi /scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta -bed /scratch/nmr15102/fhet_genome/1kb_win.bed >/scratch/nmr15102/genome_stats/basecomp/ACGT_1kb.bed