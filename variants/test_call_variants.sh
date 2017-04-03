#!/bin/bash
#SBATCH --ntasks 5
#SBATCH --cpus-per-task 4
#SBATCH -N 1-1

parallel -a /scratch/nmr15102/fhet_genome/chunks/gensplit_300.00001.bed sh ~/het_grand/variants/call_variants_subregion.sh $IN