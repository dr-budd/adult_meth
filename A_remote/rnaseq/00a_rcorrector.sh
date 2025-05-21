#!/bin/bash

# specify job checkpoint interval (-c) but only when the server executing the job is shutdown (s).
#PBS -c s
# merge standard output and standard error streams into the named file.
#PBS -j oe
# set the name of the job
#PBS -N Rcorrector
# send mail at batch job abort/exit to the email address provided.
#PBS -m ae
# advise the scheduler about the amount of physical memory required.
#PBS -l pmem=150gb
# Advise the scheduler that the job requires # cpu (ppn) from # node.
#PBS -l nodes=1:ppn=40
#Advise the scheduler that this job will have completed within 99 hours.
#PBS -l walltime=99:00:00
#Send mail at batch job abort/exit to the Email address provided.
#PBS -M jc303676@my.jcu.edu.au

#script to run Rcorrector in batch
source /etc/profile.d/modules.sh
module load jellyfish

for file in /rdsi/vol08/homes/jc303676/AGRF_RNASEQ_CAGRF18986/AGRF_CAGRF18986_CD27YANXX/RawData/R*R1.fastq;do
        R1=$file
        R2=${file/_R1/_R2}
        R1m=$(echo $file | sed 's|X_.*||')
	#echo $R1
	#echo $R2
	#echo $R1m
        nohup perl /rdsi/vol08/homes/jc303676/Software/rcorrector/run_rcorrector.pl -verbose -k 31 -t 1 -1 $R1 -2 $R2 -od /rdsi/vol08/homes/jc303676/AGRF_RNASEQ_CAGRF18986/AGRF_CAGRF18986_CD27YANXX/RawData/CD27YANXX_CorrectedReads &> ''$R1m'_Rcorrector.log'

done

#end script#

