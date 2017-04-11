#!/bin/bash
#SBATCH --ntasks 20
#SBATCH --cpus-per-task 1


# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load r/3.2.3 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

date

$parallel -a /scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.fai.10kbmin sh ~/het_grand/sweepfinder/write_sweepfinder_scaffold.sh

date
