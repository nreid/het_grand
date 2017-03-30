#!/bin/bash
#SBATCH --ntasks 10
#SBATCH --cpus-per-task 2

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load samtools 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/het_align.sh ::: $(cat /scratch/nmr15102/heteroclitus_bam/BP.list)

