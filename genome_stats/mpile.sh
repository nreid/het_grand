#!/bin/bash

module load samtools/1.3

GEN=/scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta
BAM=/scratch/nmr15102/variants/meta/bams.list

samtools mpileup -r $1 --ff 3852 --rf 2 -q 30 -Q 20 -f $GEN -b $BAM 