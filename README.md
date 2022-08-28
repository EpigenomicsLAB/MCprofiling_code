## MC profiling: a novel approach to analyze DNA methylation heterogeneity in genome-wide bisulfite sequencing data

This directory contains the R software that implements MC profiling.
Methylation Class (MC) profiling is a genome-wide approach to the study of DNA methylation heterogeneity from bulk bisulfite sequencing experiments. 
This approach is built on the concept of MCs, groups of DNA molecules sharing the same number of methylated cytosines. 
The relative abundances of MCs from sequencing reads incorporates the information on the average methylation, and directly informs on the methylation level of each molecule. 

For further details on MC profiling, see https://doi.org/10.1101/2022.07.06.498979 or contact giulia.deriso@unina.it

## Usage
## 1) Preliminary step
MC profiling is built on top of the output of the EpiStatProfiler package (https://github.com/BioinfoUninaScala/epistats).
The input of EpiStatProfiler consist in an alignment file (provided as BAM file) and a reference genome (provided as FASTA file).
Among the output files of EpiStatProfiling, the compressed matrix containing the epiallele composition for each analysed genomic region is needed to run the MC profiling code.
This epiallele matrix can be generated by running the runEpistat scripts:
```
RScripts run_Epistat/0_runEpistat.R -b <path to bam file> -g <path to reference genome fasta> -t 50 -m "CG" --modeAnalysis "n_cg" -n 4 --minBinsize 8 --maxBinsize 100 --cores <number of cores> 
```
As a first step, EpiStatProfiler will retrieve all regions holding 4 CpGs covered by at least 50 reads in the alignment file;
then, for each region it will analyze the arrangments of methylated cytosines in the sequencing reads. The arrangments found will be saved in the epiallele matrix file, along with the number of supporting reads and the coordinates of the region.
The epiallele matrix file will look like this:
```
|0111|2|chr10_294499_294532|*|
|1101|1|chr10_294499_294532*|
|1110|2|chr10_294499_294532|*|
|1111|54|chr10_294499_294532|*|
|1011|1|chr10_294510_294561|*|
|1101|2|chr10_294510_294561|*|
|1111|56|chr10_294510_294561|*|
```
