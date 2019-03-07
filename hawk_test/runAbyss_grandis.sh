#!/bin/bash

abyssDir=~/bin/abyss/bin
suffix=25_49
cutoff=49

file=case_kmers.fasta
$abyssDir/ABYSS -k25 -c0 -e0 $file -o case_abyss.fa
awk '!/^>/ { next } { getline seq } length(seq) >= 49 { print $0 "\n" seq }' case_abyss.fa > case_abyss.$suffix.fa

file=control_kmers.fasta
$abyssDir/ABYSS -k25 -c0 -e0 $file -o control_abyss.fa
awk '!/^>/ { next } { getline seq } length(seq) >= 49 { print $0 "\n" seq }' control_abyss.fa > control_abyss.$suffix.fa

