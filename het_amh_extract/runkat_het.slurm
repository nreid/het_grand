#!/bin/bash
#SBATCH --ntasks 40
#SBATCH --cpus-per-task 2

source /etc/profile.d/modules.sh

module load samtools/1.3

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/het_amh_extract/runkat_het.sh ::: $(cat /scratch/nmr15102/variants/meta/het_males.list)
$parallel sh ~/het_grand/het_amh_extract/runkat_het.sh ::: $(cat /scratch/nmr15102/variants/meta/het_females.list)
$parallel sh ~/het_grand/het_amh_extract/runkat_het.sh ::: $(cat /scratch/nmr15102/variants/meta/het_rna_females.list)
$parallel sh ~/het_grand/het_amh_extract/runkat_het.sh ::: $(cat /scratch/nmr15102/variants/meta/het_rna_males.list)
