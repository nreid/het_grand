#!/bin/bash

#1 fastq1
#2 fastq2
#3 ID table

# need to load parallel and samtools

outfile=$(echo $1 | sed 's/.sam//')

# execute bwa command line, pipe to samblaster to mark duplicates and create files containing discordant and split alignments, then to samtools to sort output. 
samtools view -S -h -u $1 | samtools sort -T $outfile - >$outfile.bam

samtools index $outfile.bam