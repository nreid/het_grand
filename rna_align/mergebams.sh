#!/bin/bash

basename=$1

inbam=/scratch/nmr15102/tolerance_rnaseq/tolerance_rnaseq
outdir=/scratch/nmr15102/tolerance_rnaseq/alignments

BAM=/home/nmr15102/bin/bamtools/bin/bamtools

find $inbam -name "*$basename*bam" >$basename.list

$BAM merge -list $basename.list -out $outdir/$basename.bam

rm $basename.list