#!/bin/bash

# 4.1 Make a text file using Stringtie run file but only running “outname” command and then copy-pasting the output into a text file

for file in /rdsi/vol08/homes/jc303676/AGRF_RNASEQ_AGRFBOTH/*.bam; do
        outname=$(echo $file | sed 's|.bam|_STRG.gtf|')
	echo $outname
done > 04.1_merge_list.txt

# 4.2 run the merge command (no need to submit a job here)

module load stringtie

nohup stringtie --merge -p 10 -G /rdsi/vol08/homes/jc303676/AGRF_RNASEQ_AGRFBOTH/TemasekBarraGenome_Annotation.gtf -o stringtie_merged.gtf /rdsi/vol08/homes/jc303676/AGRF_RNASEQ_AGRFBOTH/04.1_merge_list.txt &> 04_stringtie_merge.log &

## end script