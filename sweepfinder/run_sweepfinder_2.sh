#!/bin/bash

SF2=~/bin/SF2/SweepFinder2

indir=/scratch/nmr15102/sweepfinder

specfile=$(echo $1 | sed 's/freq/spect/')
gridfile=$(echo $1 | sed 's/\..*/.grid/')
outfile=$(echo $1 | sed 's/freq/out/')

echo $specfile $gridfile $outfile $1
$SF2 -f $indir/$1 $indir/$specfile
$SF2 -lu $indir/$gridfile $indir/$1 $indir/$specfile $indir/$outfile