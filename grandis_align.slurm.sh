#!/bin/bash
#SBATCH --ntasks 12
#SBATCH --cpus-per-task 1 

module load samtools 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

parallel sh ~/het_grand/grandis_align.sh ::: $(ls /scratch/nmr15102/gradis_seq/combined1P*) :::+  $(ls /scratch/nmr15102/gradis_seq/combined2P*) ::: /home/nmr15102/het_grand/sexes.txt
