#!/bin/bash
#SBATCH --ntasks 20
#SBATCH --cpus-per-task 2

source /etc/profile.d/modules.sh

module load samtools/1.3

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/rna_align/mergebams_osmotic.sh ::: $(find /scratch/nmr15102/rnaseq_fastas/osmotic_rnaseq/ -name "*pt.bam" | grep -oP "(?<=Sample_)[^/]+" | sed 's/AWJRDD00._//' | sort | uniq)
