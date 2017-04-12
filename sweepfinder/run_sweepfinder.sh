#!/bin/bash

SF2=~/bin/SF2/SweepFinder2

indir=/scratch/nmr15102/sweepfinder/input
griddir=/scratch/nmr15102/sweepfinder/gridfiles
outdir=/scratch/nmr15102/sweepfinder/output
specdir=/scratch/nmr15102/sweepfinder/globalspec

specfile=$(echo $1 | grep -oP "(?<=\.1_)[A-Z]+_[a-z]+").spect
gridfile=$(echo $1 | grep -oP ".*.1").grid
outfile=$(echo $1 | sed 's/geno/out/')

$echo $specfile $gridfile $outfile $1
$SF2 -lu $griddir/$gridfile $indir/$1 $specdir/$specfile $outdir/$outfile