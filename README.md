## MC profiling: a novel approach to analyze DNA methylation heterogeneity in genome-wide bisulfite sequencing data

This directory contains the R software that implements MC profiling.
Methylation Class (MC) profiling is a genome-wide approach to the study of DNA methylation heterogeneity from bulk bisulfite sequencing experiments. 
This approach is built on the concept of MCs, groups of DNA molecules sharing the same number of methylated cytosines. 
The relative abundances of MCs from sequencing reads incorporates the information on the average methylation, and directly informs on the methylation level of each molecule. 

For further details on MC profiling, see https://doi.org/10.1101/2022.07.06.498979 or contact giulia.deriso@unina.it

## Usage
## 1) Preliminary step
In order to perform MC profiling, we need to run the EpiStatProfiler package (https://github.com/BioinfoUninaScala/epistats).
```
RScripts 0_runEpistat.R -b <path to bam file> -g <path to reference genome fasta> -t 50 -m "CG" --modeAnalysis "n_cg" -n 4 --minBinsize 8 --maxBinsize 100 --cores <number of cores> 
```
EpiStatProfiler takes in input an alignment file (provided as BAM file) and a reference genome (provided as FASTA file).
As a first step, EpiStatProfiler will retrieve all regions holding 4 CpGs (epiloci) covered by at least 50 reads in the alignment file;
then, for each region it will analyze the arrangments of methylated cytosines in the sequencing reads. The arrangments found will be saved in the epiallele matrix file, along with the number of supporting reads and the coordinates of the region.
The epiallele matrix file ("<sample>"_epiAnalysis.txt) will look like this:

| var1 | Freq | id | strand |
|:----:|:-:|:-------------------:|:-:|
| 0111 | 2 | chr10_294499_294532 | * |
| 1101 |1  | chr10_294499_294532 | * |
| 1110 |2  | chr10_294499_294532 | * |
| 1111 |54 | chr10_294499_294532 | * |
| 1011 | 1 | chr10_294510_294561 | * |
| 1101 | 2 | chr10_294510_294561 | * |
| 1111 | 56| chr10_294510_294561 | * |

EpiStatProfiler also compute some statistics for the analzed regions, which will be saved in a bed file.
The bed file ("<sample>"_intervals.bed) will look like this:

| seqnames | start	| end	| width	| strand	| dist	| epi	| singleton	| maxfreq	| shannon	| mean_met	| num_cg |	num_reads |
|:--------:|:------:|:---:|:-----:|:-------:|:-----:|:---:|:---------:|:-------:|:-------:|:---------:|:------:|:----------:|
| chr10	| 294499	| 294532	| 34	| *	| 2.67	| 4	| 1	| 1111	| 0.55	| 0.98	| 4	| 59 |
| chr10	| 294510	| 294561	| 52	| *	| 4.17	| 3	| 1	| 1111	| 0.34	| 0.99	| 4	| 59 |

For further details on the EpiStatProfiler package, see https://github.com/BioinfoUninaScala/epistats/

## 2) Filtering of epiloci
By executing the script I_epilociFiltering.R, we will select the epiloci with coverage higher or equal to 50 reads and lower than 99th percentile. In case of overlap, we will retain only one of the epiloci.
For each sample, two new files will be generated: "<sample>"_epiAnalysis_filtered.txt and "<sample>"_intervals_filtered.bed.
Note that the script can deal with multiple samples, provided that the epiallele matrix and the bed file are in the same input directory and that they can be matched through their "<sample>" prefix. The input dir can be set by the user changing the value of the inputDir variable in the script.

## 3) Computation of MC profiles

  


