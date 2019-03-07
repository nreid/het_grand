#!/bin/bash
#SBATCH --ntasks 40
#SBATCH --cpus-per-task 2

source /etc/profile.d/modules.sh

module load samtools/1.3

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

$parallel sh ~/het_grand/grandis_extract/runkat_grandis.sh ::: $(cat /scratch/nmr15102/variants/meta/grand_males.list)
$parallel sh ~/het_grand/grandis_extract/runkat_grandis.sh ::: $(cat /scratch/nmr15102/variants/meta/grand_females.list)
$parallel sh ~/het_grand/grandis_extract/runkat_grandis.sh ::: $(cat /scratch/nmr15102/variants/meta/grand_rnaseq.list)
