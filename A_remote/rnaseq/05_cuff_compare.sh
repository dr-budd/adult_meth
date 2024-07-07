#!/bin/bash

# Examine how the transcripts compare with the reference annotation (optional)
# NB: using cuffcompare because gff utilities was not installed on the HPC

Module load cufflinks/2.2.1

cuffcompare –r TemasekBarraGenome_Annotation.gtf stringtie_merged.gtf

## end script


