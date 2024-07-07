#!/bin/bash

## run fastqc for all WGBS files

echo "working directory: "$PWD

sftloc=$PWD/../../Software/FastQC/fastqc	## software location
datadir=$PWD/../rawdata
outdir=$PWD/../analysis/fastqc

nohup $sftloc $datadir/*fastq.gz --outdir=$outdir &>fastqc_nohup.out

## end script
