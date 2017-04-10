#!/bin/bash
#SBATCH --ntasks 40
#SBATCH --cpus-per-task 1


# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load r/3.2.3 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

date
$parallel -a /scratch/nmr15102/fhet_genome/chunks/gensplit_300.00009.bed sh ~/het_grand/variants/haploidize_subregion.sh
date
