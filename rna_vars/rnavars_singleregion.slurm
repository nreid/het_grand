#!/bin/bash
#SBATCH --ntasks 100
#SBATCH --cpus-per-task 1

source /etc/profile.d/modules.sh

module load samtools/1.3
module load xz

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

BED=/scratch/nmr15102/tolerance_rnaseq/callregions_rna_intersect.bed
$parallel bash ~/het_grand/rna_vars/rnavars_singleregion.sh ::: $(cat $BED | cut -f 1-4 | sed 's/\t/:/' | sed 's/\t/-/')
