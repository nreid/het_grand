#!/bin/bash


BGZ=~/bin/htslib/bgzip

sh ~/het_grand/variants/merge_subvcfs.sh | $BGZ -c >/scratch/nmr15102/variants/hetgrand.vcf.gz