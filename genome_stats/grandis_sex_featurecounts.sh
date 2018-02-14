#!/bin/bash

FC=~/bin/subread-1.5.3-Linux-x86_64/bin/featureCounts
ANN=/scratch/nmr15102/fhet_genome/1kbwin.saf
OUT=/scratch/nmr15102/genome_stats/grandsex/1kb_sex_counts.txt
LIST=/scratch/nmr15102/variants/meta/bams.list

$FC \
-F SAF \
-Q 30 \
--primary \
-O \
-p \
-B \
-a $ANN -o $OUT $(cat $LIST | tr "\n" " ")
