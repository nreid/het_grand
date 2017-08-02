#!/bin/bash
#SBATCH --ntasks 10
#SBATCH --cpus-per-task 2
#SBATCH -N 1-1


source /etc/profile.d/modules.sh

module load samtools/1.3

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/rna_align/rna_align.sh ::: $(find /scratch/nmr15102/tolerance_rnaseq/ -name "*fastq.gz")
