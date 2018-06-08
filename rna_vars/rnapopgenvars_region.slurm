#!/bin/bash
#SBATCH --ntasks 100
#SBATCH --cpus-per-task 1

source /etc/profile.d/modules.sh

module load samtools/1.3
module load xz

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/rna_vars/rnapopgenvars_region.sh ::: $(cat /scratch/nmr15102/tolerance_rnaseq/callregions_rna_union.bed | cut -f 1 | uniq)
