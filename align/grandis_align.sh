#!/bin/bash

#1 fastq1
#2 fastq2
#3 ID table

# need to load parallel and samtools


sam=$(echo $1 | grep -oP "BU[0-9]+")
pop=$(cat $3 | grep $sam | cut -f 3)

fq1=$1
fq2=$2

rg=$(echo \@RG\\tID:$sam.combined\\tPL:Illumina\\tPU:x\\tLB:combined\\tSM:$sam.$pop)

outdir=/scratch/nmr15102/grandis_align
outroot=$sam\_$pop



bwagenind=/scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.bgz

BWA=~/bin/bwa/bwa
SBL=~/bin/samblaster/samblaster

cmdline=$BWA\ mem\ $bwagenind\ -t\ 2\ -R\ $rg\ $fq1\ $fq2
echo $cmdline

# execute bwa command line, pipe to samblaster to mark duplicates and create files containing discordant and split alignments, then to samtools to sort output. 
$cmdline | $SBL -e -d $outdir/$outroot.disc.sam -s $outdir/$outroot.split.sam | samtools view -S -h -u - | samtools sort -T $outdir/$outroot - >$outdir/$outroot.bam

