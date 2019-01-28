#!/bin/bash
#SBATCH --ntasks 40
#SBATCH --cpus-per-task 2


source /etc/profile.d/modules.sh

module load samtools/1.3

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/rna_align/rna_align_grandis.sh ::: $(find /scratch/nmr15102/grandis_rnaseq/fastqs -name "*fq.gz" | grep R1)
