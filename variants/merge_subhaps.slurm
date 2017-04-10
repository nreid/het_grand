#!/bin/bash


BGZ=~/bin/htslib/bgzip

sh ~/het_grand/variants/merge_subhaps.sh | $BGZ -c >/scratch/nmr15102/variants/hetgrand.hap.vcf.gz