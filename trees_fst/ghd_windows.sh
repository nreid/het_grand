#!/bin/bash

source /etc/profile.d/modules.sh 

module load r/3.3.3

TAB=~/bin/htslib/tabix

VCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz

reg=$(echo $1 | sed 's/:/	/' | sed 's/-/	/')
$TAB $VCF $1 | Rscript ~/het_grand/trees_fst/ghd_stream.R | sed "s/^/$reg	/" >/scratch/nmr15102/intdist/singles/$1.ghd.out