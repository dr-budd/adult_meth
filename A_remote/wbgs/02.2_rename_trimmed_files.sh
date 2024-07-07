#!/bin/bash

## navigate to trimming dir
cd $PWD/../analysis/trimming

echo "location: "$PWD

## format
## find . -type f -name '*YOUR_IDENTIFIER_HERE*' | while read FILE ; do
##    newfile="$(echo ${FILE} |sed -e 's/OLDBIT/NEWBIT/')" ;
##    mv "${FILE}" "${newfile}" ;
## done 

## R1s
find . -type f -name '*_1.fq.gz' | while read FILE; do 
	newfile="$(echo ${FILE} |sed -e 's/val_1/trimmed/')" ;
    mv "${FILE}" "${newfile}" ;
done

## R2s
find . -type f -name '*_2.fq.gz' | while read FILE; do
        newfile="$(echo ${FILE} |sed -e 's/val_2/trimmed/')" ;
    mv "${FILE}" "${newfile}" ;
done

## end script
