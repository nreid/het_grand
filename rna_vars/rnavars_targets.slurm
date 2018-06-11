#!/bin/bash
#SBATCH --ntasks 100
#SBATCH --cpus-per-task 1

source /etc/profile.d/modules.sh

module load samtools/1.3
module load xz

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/rna_vars/rnavars_targets.sh ::: $(seq 1 50 194207)
