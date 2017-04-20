#!/bin/bash

# "module" alone was an unrecognized command. very confusing, the parallel slurm script runs this before "module" though. so I'm trying it. 
source /etc/profile.d/modules.sh

module load r/3.2.3 


TAB=~/bin/htslib/tabix
BGZ=~/bin/htslib/bgzip

zcat /scratch/nmr15102/variants/hetgrand.hap.vcf.gz | \
grep -v ^# |  \
Rscript ~/het_grand/trees_fst/fst_stream_bysite.R | \
$BGZ -c >het_grand_fst.gz