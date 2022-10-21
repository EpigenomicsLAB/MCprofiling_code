## MC profiling: a novel approach to analyze DNA methylation heterogeneity in genome-wide bisulfite sequencing data

This directory contains the R software that implements MC profiling.
Methylation Class (MC) profiling is a genome-wide approach to the study of DNA methylation heterogeneity from bulk bisulfite sequencing experiments. 
This approach is built on the concept of MCs, groups of DNA molecules sharing the same number of methylated cytosines. 
The relative abundances of MCs from sequencing reads incorporates the information on the average methylation, and directly informs on the methylation level of each molecule. 

For further details on MC profiling, see https://doi.org/10.1101/2022.07.06.498979 or contact giulia.deriso@unina.it

## Set up R environment
You can set up the R environment (R version >= 3.6) to run the MC profiling code through anaconda (https://www.anaconda.com/products/distribution), following this steps:
```
conda create -n MCprofiling r-base r-essentials
conda activate MCprofiling
conda install -c r r-tidyverse
git clone https://github.com/EpigenomicsLAB/MCprofiling_code.git
cd MCprofiling_code
R
source("install_dependencies.R")
```
## Usage
## 1) Preliminary step
In order to perform MC profiling, we need to run the EpiStatProfiler package (https://github.com/BioinfoUninaScala/epistats).
```
Rscript 0_runEpistat.R -b <path to bam file> -g <path to reference genome fasta> -t 50 -m "CG" --modeAnalysis "n_cg" -n 4 --minBinsize 8 --maxBinsize 100 --cores <number of cores> 
```
EpiStatProfiler takes in input an alignment file (provided as BAM file) and a reference genome (provided as FASTA file).
As a first step, EpiStatProfiler will retrieve all regions holding 4 CpGs (epiloci) covered by at least 50 reads in the alignment file;
then, for each region it will analyze the arrangments of methylated cytosines in the sequencing reads. The arrangments found will be saved in the epiallele matrix file, along with the number of supporting reads and the coordinates of the region.
The epiallele matrix file (\<sample\>_epiAnalysis.txt) will look like this:

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
The bed file (\<sample\>_intervals.bed) will look like this:

| seqnames | start	| end	| width	| strand	| dist	| epi	| singleton	| maxfreq	| shannon	| mean_met	| num_cg |	num_reads |
|:--------:|:------:|:---:|:-----:|:-------:|:-----:|:---:|:---------:|:-------:|:-------:|:---------:|:------:|:----------:|
| chr10	| 294499	| 294532	| 34	| *	| 2.67	| 4	| 1	| 1111	| 0.55	| 0.98	| 4	| 59 |
| chr10	| 294510	| 294561	| 52	| *	| 4.17	| 3	| 1	| 1111	| 0.34	| 0.99	| 4	| 59 |

For further details on the EpiStatProfiler package, see https://github.com/BioinfoUninaScala/epistats/

## 2) Filtering of epiloci
By executing the script I_epilociFiltering.R, we will select the epiloci with coverage higher or equal to 50 reads and lower than 99th percentile. In case of overlap, we will retain only one of the epiloci.
For each sample, two new files will be generated: \<sample\>_epiAnalysis_filtered.txt and "<sample>"_intervals_filtered.bed.
Note that the script can deal with multiple samples, provided that the epiallele matrix and the bed file are in the same input directory and that they can be matched through their \<sample\> prefix. The input directory can be set by the user changing the value of the inputDir variable in the script.

## 3) Computation of MC profiles
By executing the script II_MCprofiling.R, we will compute the MC profiles for the selected epiloci.
The epialleles with the same number of methylated cytosines (indicated as 1) will be grouped in Methylation Classes, and the number of supporting reads will be summed. Then, for each possible MC, the relative abundance (i.e. the number of supporting reads out of the number of reads mapped to the epilocus) will be computed.
The MC profiles will be saved in a file named \<sample\>_MCprofiles.txt.
| id	| 0	| 1	| 2 |	3	| 4 |
|:---:|:-:|:-:|:-:|:-:|:-:|
| chr10_1000256_1000305	| 0	| 0	| 0	| 0.08064516129032258	| 0.9193548387096774 |
| chr10_100227667_100227709 |	1	| 0	| 0	| 0	| 0 |

Note that the script can deal with multiple samples if multiple epiallele filtered matrix files are provided in the same input directory. Also in this case, the input directory can be set by the user changing the value of the inputDir variable in the script.

## 4) Comparison of MC profiles
To quantify the dissimilarity between two MC profiles, we adopted the Jensen-Shannon distance (JSD).
```
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log2(x/y), na.rm = T)
```
where x and y are the vectors of the abundances of the two MC profiles to be compared.

## 5) Multiple samples handling
When multiple samples are available for a given condition, we can compute an average profile by running the script III_sampleConsensus.R.
First, we will select the epiloci shared by all the samples.
Then, for each epilocus, we will compare the MC profile among all the possible sample pairs. Only the stable epiloci (i.e., epiloci with JSD below the established threshold of 0.26) will be selected.
Finally, for each selected epilocus, we will compute the consensus MC profile by averaging the relative abundances of MCs in the samples.
The consensus profiles, along with the JSD observed in the sample pairs, will be saved in the consensus_MCprofiles.txt file.

| id	| sample1_sample2	| sample1_sample3	| sample2_sample3	| 0	| 1	| 2	| 3	| 4 |
|:---:|:---------------:|:---------------:|:---------------:|:-:|:-:|:-:|:-:|:-:|
| chr10_100227667_100227709	| 0.19892611694951864	| 0.17512680621292145	| 0.02849152161736604	| 0.9543589743589744	| 0.045641025641025644	| 0	| 0	| 0 |
| chr10_100992192_100992223	| 0.14617785097801772	| 0.19261120760552786	| 0.1715471175232086	| 0.8390612202268067 |	0.12982019272525383	| 0.011494252873563216	| 0.013877207737594613	| 0.0057471264367816065 |
  
Note that the script will consider all the samples for which an MCprofiles.txt file is provided in the input directory.
## 6) MC profile classification
To further improve the interpretability of MC profiles, we can assigne each MC profile to a Methylation Patterns by running the script IV_MPClassification.R (see https://doi.org/10.1101/2022.07.06.498979 for details).
Each MC profile is compared to 5 archetypal profiles by computing the JSD. The profile is then assigned to the MP which correspond to the archetype with the lowest JSD.
The MP, along with the JSD from the 5 archetypes, is saved in the  MP_classification file.
Here is an example of classification for the previously computed consensus profiles.

| id	| sample1_sample2	| sample1_sample3	| sample2_sample3	| 0	| 1	| 2	| 3	| 4	| D1	| D2	| D3	| D4	| D5	| Class |
|:---:|:---------------:|:---------------:|:---------------:|:-:|:-:|:-:|:-:|:-:|:---:|:---:|:---:|:---:|:---:|:-----:|
| chr10_100227667_100227709	| 0.19892611694951864	| 0.17512680621292143	| 0.02849152161736604	| 0.9543589743589744	| 0.04564102564102564	| 0 |	0	| 0	| 0.2065395137803402	| 1	| 0.5846289702431131	| 0.9530208819678299	| 0.728732604086668	| D1 |
| chr10_100992192_100992223	| 0.14617785097801772	| 0.19261120760552783	| 0.1715471175232086	| 0.8390612202268067	| 0.12982019272525383	| 0.01149425287356321	| 0.01387720773759461	| 0.0057471264367816 |	0.146021950545462	| 0.9686685813177589	| 0.521633692547548	| 0.863221494527554	| 0.620654534924665	| D1 |
