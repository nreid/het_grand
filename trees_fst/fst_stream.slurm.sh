#!/bin/bash
#SBATCH --ntasks 60
#SBATCH --cpus-per-task 1


# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load r/3.2.3 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

INTERVAL=/scratch/nmr15102/sweepfinder/null_invervals_10kb.txt
DIR=/scratch/nmr15102/sweepfinder
OUT=fst.null.out

date

$parallel -a $INTERVAL sh ~/het_grand/trees_fst/fst_stream.sh >$DIR/$OUT

date
