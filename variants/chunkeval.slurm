#!/bin/bash
#SBATCH --ntasks 40
#SBATCH --cpus-per-task 1



source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/variants/chunkeval.sh ::: $(ls /scratch/nmr15102/variants/subvars/*gz | head -n 300) >chunkeval_out.txt
