# :fish::dna: Effects of temperature on DNA methylation and sex in adult barramundi

## :page_facing_up: Associated Publication

[![Journal](https://img.shields.io/badge/Published_in-Aquaculture-6E1E78.svg)](https://doi.org/10.1016/j.aquaculture.2025.742882)

> **Budd, A.M.**, Huerlimann, R., Guppy, J.L., Pinto, R.C.C., Domingos, J.A. & Jerry, D.R. (2026).  
> Mechanisms driving temperature-induced early sex change in barramundi (*Lates calcarifer*).  
> *Aquaculture*, **610**, 742882.  
> [https://doi.org/10.1016/j.aquaculture.2025.742882](https://doi.org/10.1016/j.aquaculture.2025.742882)

---

## :file_folder: Structure

The code is separated into the following two folders:

`A_remote` - Contains all scripts run remotely. Note: These scripts were designed for a HPC facility using a PBS batch-queue system. Run times and memory usage allocations may need to be adjusted for other platforms.

`B_local` - Contains all scripts run locally. 

## :chart_with_upwards_trend: Data

Raw sequence data and metadata are available in the NCBI Sequence Read Archive under BioProject accession PRJNA1263116, with BioSample accessions SAMN48511396-SAMN48511415. 
SRA accessions for RNA-seq data are SRR33638345-SRR33638356, and for WGBS data are SRR33696300 - SRR33696311.

## :woman_technologist: Author
Alyssa Budd (alyssa.budd@my.jcu.edu.au)

## :bouquet: Acknowledgements
Scripts were developed in collaboration with Queensland Facility for Advanced Bioinformatics (QFAB). Specifically, Mike Thang wrote the original code for WGBS methlylation calling and Anne Bernard drafted the code for running multiple comparisons between treatments for the RNAseq count data.

## :copyright: License
This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE.txt) file for details.