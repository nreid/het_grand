#!/bin/bash
#SBATCH --ntasks 5
#SBATCH --cpus-per-task 4
#SBATCH -N 1-1


# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load samtools 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

date
$parallel -a /scratch/nmr15102/fhet_genome/chunks/gensplit_300.00062.bed sh ~/het_grand/variants/call_variants_subregion.sh
date
