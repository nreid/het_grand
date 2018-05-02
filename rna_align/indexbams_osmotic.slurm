#!/bin/bash
#SBATCH --ntasks 40
#SBATCH --cpus-per-task 1

source /etc/profile.d/modules.sh

module load samtools/1.3

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/rna_align/indexbams.sh ::: $(find /scratch/nmr15102/osmotic_rnaseq/alignments -name "*bam")
