#!/bin/bash

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh 

module load samtools/1.3

for file in $(cat /scratch/nmr15102/variants/meta/bams.list)
do samtools stats $file | gzip -c >$file.stats
done