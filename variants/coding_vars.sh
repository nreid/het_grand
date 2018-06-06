#!/bin/bash

source /etc/profile.d/modules.sh
module load xz

# script gets and annotates only variants residing in CDS

#programs
BED=/home/nmr15102/bin/bedtools2/bin/bedtools
VAP=/home/nmr15102/bin/vcflib/bin/vcfallelicprimitives
BGZ=/home/nmr15102/bin/htslib/bgzip
EFF=/home/nmr15102/bin/snpEff/snpEff.jar

#files and directories
GFF=/home/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.gff
VCF=/scratch/nmr15102/variants/hetgrand.hap.vcf.gz
OUTDIR=/scratch/nmr15102/variants


 $BED intersect -u -header -a $VCF -b <(awk '$3 ~ /CDS/' $GFF) | \
 $VAP | \
 $BED intersect -u -header -a stdin -b <(awk '$3 ~ /CDS/' $GFF) | \
 sed 's/:.:.:.:.:.:.:.//g' | \
 java -jar $EFF ann killifish - | \
 $BGZ >hetgrand.hap_CDSann.vcf.gz