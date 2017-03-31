#!/bin/bash
#SBATCH --ntasks 10
#SBATCH --cpus-per-task 2
#SBATCH -N 1-1

source /etc/profile.d/modules.sh

module load samtools 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/align/het_align.sh ::: $(cat /scratch/nmr15102/heteroclitus_bam/F.list)

