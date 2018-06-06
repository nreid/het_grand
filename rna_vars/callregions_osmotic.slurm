#!/bin/bash
#SBATCH --ntasks 100
#SBATCH --cpus-per-task 1

source /etc/profile.d/modules.sh

module load samtools/1.3

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/rna_vars/callregions_osmotic.sh ::: $(cat /scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.fai | cut -f 1)
