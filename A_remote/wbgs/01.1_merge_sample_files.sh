#!/bin/bash

# to make this list (very round about, sorry) 
# step 1
# ls -1 *R1*.gz | awk -F '_' '{print $1}' | sort | uniq > ID
# step 2
# for i in `cat ./ID` ; do echo cat $i\*_R1.fastq.gz ">" $i\_R1.merged.fastq.gz ;  echo cat $i\*_R2.fastq.gz ">" $i\_R2.merged.fastq.gz ; done
# step 3
# copy paste the output into here

cat 24D1*_R1.fastq.gz > 24D1_R1.merged.fastq.gz
cat 24D1*_R2.fastq.gz > 24D1_R2.merged.fastq.gz
cat 24D2*_R1.fastq.gz > 24D2_R1.merged.fastq.gz
cat 24D2*_R2.fastq.gz > 24D2_R2.merged.fastq.gz
cat 24D3*_R1.fastq.gz > 24D3_R1.merged.fastq.gz
cat 24D3*_R2.fastq.gz > 24D3_R2.merged.fastq.gz
cat 28D1*_R1.fastq.gz > 28D1_R1.merged.fastq.gz
cat 28D1*_R2.fastq.gz > 28D1_R2.merged.fastq.gz
cat 28D2*_R1.fastq.gz > 28D2_R1.merged.fastq.gz
cat 28D2*_R2.fastq.gz > 28D2_R2.merged.fastq.gz
cat 28D3*_R1.fastq.gz > 28D3_R1.merged.fastq.gz
cat 28D3*_R2.fastq.gz > 28D3_R2.merged.fastq.gz
cat 34D1*_R1.fastq.gz > 34D1_R1.merged.fastq.gz
cat 34D1*_R2.fastq.gz > 34D1_R2.merged.fastq.gz
cat 34D2*_R1.fastq.gz > 34D2_R1.merged.fastq.gz
cat 34D2*_R2.fastq.gz > 34D2_R2.merged.fastq.gz
cat 34D5*_R1.fastq.gz > 34D5_R1.merged.fastq.gz
cat 34D5*_R2.fastq.gz > 34D5_R2.merged.fastq.gz
cat FD1*_R1.fastq.gz > FD1_R1.merged.fastq.gz
cat FD1*_R2.fastq.gz > FD1_R2.merged.fastq.gz
cat FD2*_R1.fastq.gz > FD2_R1.merged.fastq.gz
cat FD2*_R2.fastq.gz > FD2_R2.merged.fastq.gz
cat FD3*_R1.fastq.gz > FD3_R1.merged.fastq.gz
cat FD3*_R2.fastq.gz > FD3_R2.merged.fastq.gz

# to run:
# chmod +x this_files_name.sh
# nohup ./this_files_name.sh &

# end script
