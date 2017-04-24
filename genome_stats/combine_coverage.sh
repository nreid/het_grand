#!/bin/bash

BED=~/bin/bedtools2/bin/bedtools
BGZ=~/bin/htslib/bgzip

zcat admix.coverage.gz | awk '$3>9' | awk '{OFS="\t"}{start=$2-1}{print $1,start,$2,$3}'

$BED intersect -sorted -a <(zcat admix.coverage.gz | awk '$3>9' | awk '{OFS="\t"}{start=$2-1}{print $1,start,$2,$3}') -b <(zcat north.coverage.gz | awk '$3>9' | awk '{OFS="\t"}{start=$2-1}{print $1,start,$2,$3}') | \
$BED intersect -sorted -a stdin -b <(zcat south.coverage.gz | awk '$3>9' | awk '{OFS="\t"}{start=$2-1}{print $1,start,$2,$3}') | \
$BED intersect -sorted -a stdin -b <(zcat grand.coverage.gz | awk '$3>19' | awk '{OFS="\t"}{start=$2-1}{print $1,start,$2,$3}') | \
$BGZ -c | cat >coverage.combined.gz