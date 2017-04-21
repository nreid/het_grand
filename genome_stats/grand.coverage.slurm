#!/bin/bash

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load samtools/1.3

LIST=/scratch/nmr15102/variants/meta/grand.list

sh ~/het_grand/genome_stats/coverage_bysite.sh $LIST