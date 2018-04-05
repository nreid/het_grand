#!/bin/bash

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh 

source ~/bin/parallel-slurm/parallel-slurm-setup.sh

module load r/3.3.3

TAB=~/bin/htslib/tabix

VCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz

for num in {1..18873}
do win=$(sed -n $(echo $num)p /home/nmr15102/fhet_genome/100kb_win.bed)
reg=$(echo $win | sed 's/ /:/' | sed 's/ /-/')
tre=$(echo $win | sed 's/ /	/g')
$TAB $VCF $reg | Rscript ~/het_grand/trees_fst/ghd_stream.R | sed "s/^/$tre	/"
done