#!/bin/bash

#where to find module load command
source /etc/profile.d/modules.sh

#load hisat2 and samtools
module load hisat2
module load samtools/1.5

for file in $PWD/*.sam; do
        samfile=$file
        outname=$(echo $file | sed 's|.sam|.bam|')

#echo $samfile
#echo $outname

nohup samtools sort -@ 10 -o $outname $samfile &> /rdsi/vol08/homes/jc303676/AGRF_RNASEQ_AGRFBOTH/02_samtools_view_sort.log &

done

# end script #