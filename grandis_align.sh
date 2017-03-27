#!/bin/bash

#1 fastq1
#2 fastq2
#3 ID table


sam=$(echo $1 | grep -oP "BU[0-9]+")
pop=$(cat $3 | grep $sam | cut -f 2)

fq1=$1
fq2=$2

rg=$(echo \@RG\\tID:$sam.combined\\tPL:Illumina\\tPU:x\\tLB:combined\\tSM:$sam.$pop)


# module load 

# bwagenind=/home/nreid/popgen/kfish3/killifish20130322asm.fa

# BWA=~/bin/bwa/bwa
# SBL=~/bin/samblaster/samblaster

# cmdline=$BWA\ mem\ $bwagenind\ -t\ 24\ -R\ $rg\ $fq1\ $fq2
# echo $cmdline

# #execute bwa command line, pipe to samblaster to mark duplicates and create files containing discordant and split alignments, then to samtools to sort output. 
# $cmdline | $SBL -e -d $bwadir/$outroot.disc.sam -s $bwadir/$outroot.split.sam | samtools view -S -h -u - | samtools sort - $bwadir/$outroot

