#!/bin/bash
#SBATCH --ntasks 80
#SBATCH --cpus-per-task 1


# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel -a ER_in.txt sh ~/het_grand/sweepfinder/run_sweepfinder.sh
