#!/bin/bash

basename=$1

inbam=/scratch/nmr15102/rnaseq_fastas/osmotic_rnaseq
outdir=/scratch/nmr15102/osmotic_rnaseq/alignments

BAM=/home/nmr15102/bin/bamtools/bin/bamtools

find $inbam -name "*$basename*bam" >$basename.list

$BAM merge -list $basename.list -out $outdir/$basename.bam

rm $basename.list