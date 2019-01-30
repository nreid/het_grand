#!/bin/bash

FC=~/bin/subread-1.5.3-Linux-x86_64/bin/featureCounts
ANN=~/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.GeneIDs.gff
OUT=/scratch/nmr15102/grandis_rnaseq/counts/counts.txt
LIST=/scratch/nmr15102/grandis_rnaseq/alignments/meta/bams.list

$FC \
-F GFF \
-g gene \
-Q 30 \
--primary \
-p \
-a $ANN -o $OUT $(cat $LIST | tr "\n" " ")
