#!/bin/bash
#SBATCH --ntasks 12
#SBATCH --cpus-per-task 1 

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load samtools 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

parallel sh ~/het_grand/grandis_align.sh ::: $(ls /scratch/nmr15102/gradis_seq/combined1P*) :::+  $(ls /scratch/nmr15102/gradis_seq/combined2P*) ::: /home/nmr15102/het_grand/sexes.txt
