#!/bin/bash
#SBATCH --ntasks 50
#SBATCH --cpus-per-task 2

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load samtools 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel -n0 sh ~/het_grand/het_align.sh ::: $(ls /scratch/nmr15102/heteroclitus_bam/BP-*bam)

