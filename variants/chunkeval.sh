#!/bin/bash

BLAH=$(zcat $1 | \
grep -v ^# | \
cut -f 1-2 | \
awk '{OFS="\t"}{k+=1}{if (k==1) {st=$2 ; sc=$1}} END {en=$2; 	print sc,st,en}')
echo $1 $BLAH
