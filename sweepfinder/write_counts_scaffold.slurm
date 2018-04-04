#!/bin/bash
#SBATCH --ntasks 60
#SBATCH --cpus-per-task 1


# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load r/3.3.3 
module load xz

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

date

$parallel -a ~/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.fai.10kbmin sh ~/het_grand/sweepfinder/write_counts_scaffold.sh

date
