# :fish::dna: Effects of temperature on DNA methylation and sex in adult barramundi

Analysis pipeline for:

>Budd et al. (*submitted*) Mechanisms driving temperature-induced early sex change in barramundi (*Lates calcarifer*)

---

## :file_folder: Structure

The code is separated into the following two folders:

`A_remote` - Contains all scripts run remotely. Note: These scripts were designed for a HPC facility using a PBS batch-queue system. Run times and memory usage allocations may need to be adjusted for other platforms.

`B_local` - Contains all scripts run locally. 

## :chart_with_upwards_trend: Data

Raw sequence data and metadata will be archived in the Dryad data repository upon acceptance of the manuscript.

## :woman_technologist: Author
Alyssa Budd (alyssa.budd@my.jcu.edu.au)

## :bouquet: Acknowledgements
Scripts were developed in collaboration with Queensland Facility for Advanced Bioinformatics (QFAB). Specifically, Mike Thang wrote the original code for WGBS methlylation calling and Anne Bernard drafted the code for running multiple comparisons between treatments for the RNAseq count data.

## :copyright: License
This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE.txt) file for details.