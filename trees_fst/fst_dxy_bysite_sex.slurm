#!/bin/bash
#SBATCH --ntasks 80
#SBATCH --cpus-per-task 2


# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh 

module load r/3.3.3

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel -a /home/nmr15102/fhet_genome/intervals_1mb.txt sh ~/het_grand/trees_fst/fst_dxy_bysite_sex.sh
