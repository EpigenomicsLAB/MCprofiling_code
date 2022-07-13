#############################################################################################################
#  Copyright (C) 2020-21 Lab Cocozza, University of Naples Federico II, Naples
## CONFIDENTIAL
## This document is proprietary and confidential. No part of this document may be disclosed in any manner
## to a third party without the prior written consent
#############################################################################################################
#############################################################################################################

####This script filters out apiloci based on their coverage
#### and randomly selects one epilocus among overlapping epiloci.
#### In this way, only non overlapping epiloci covered by at least 50 reads will be retained for further analysis

####################################################################################
library(tidyverse)
library(parallel)
####Parameters
wd="."
setwd(wd)
inputDir="test_data"
# set number of cores
n_cores=4
# the input files should be the intervals.bed file created by EpiTool
inputFiles= list.files(inputDir) %>% .[matches("intervals.bed$", vars=.)]
sample_names= sapply(strsplit(inputFiles, split="_"), "[[",1)
inputFiles=split(inputFiles, f=sample_names)
# the input files should be the epiAnalysis.txt file created by EpiTool. If more than one file is provided, items should be placed in the same order than inputFiles
matrixFiles= list.files(inputDir) %>% .[matches("epiAnalysis.txt$", vars=.)]
sample_names= sapply(strsplit(matrixFiles, split="_"), "[[",1)
matrixFiles=split(matrixFiles, f=sample_names)

matrixFiles=matrixFiles[match(names(matrixFiles), names(inputFiles))]

# set output dir
outDir=paste0(wd, "/out_I")
#minimum coverage required for a region to be retained
min_readNum=50

#minimum coverage required for a region to be retained
#regions with coverage above the 99.9 percentiles will be fileterd out
#the user can modify the script behavior by setting a value for maxreadNum 
#the custom values have to be provided as a list of the same length of inputFiles
maxreadNum=list()

####Custom functions
removeOverlaps=function(overlap, cores)
{
  overlap=split(overlap,f = overlap$seqnames)
  overlap=map(overlap, function(x) arrange(x, start))
  overlap=mclapply(overlap, function(x) {
    i=1
    retain=c()
    while(i<=nrow(x))
    {
      retain=append(retain,i)
      chr=x[i,2]
      end=x[i,4]
      df=x[c(i:nrow(x)),]
      Next=min(which(df$seqnames==chr & df$start>end))
      i=i+Next-1
    }
    return(x[retain,])
  },mc.cores = 4)
  
  return(Reduce(bind_rows,overlap))
}

#######main
#1) read input files
intervals=map(inputFiles, function(x) read_tsv(paste(inputDir,x,sep="/")))
#2) assign region id (Chr_start_end)
intervals=map(intervals, function(x) unite(x, id, c("seqnames":"end"), remove = F, sep = "_"))
#3) filters regions with coverage below minThresh
if(length(maxreadNum)==0)
{
  maxreadNum=map(intervals, function(x) unname(quantile(x$num_reads, seq(0,1,0.001))[1000]))
}

intervals=Map(function(x,y) x %>% filter(num_reads<= y) %>% return(),
              intervals, maxreadNum)
#4) filters out regions above maxThresh
intervals=map(intervals, function(x) return(x[x$num_reads >=min_readNum,]))
#5) select non overlapping regions
intervals=map(intervals, function(x) removeOverlaps(as.data.frame(x), n_cores))
#7) load epiMatrix 
epiMatr=map(matrixFiles, function(x) read_tsv(paste(inputDir,x,sep="/"), col_types = c("c","n","c","c")))
#8) get epialleles of selected regions in intervals
epiMatr=Map(function(x,y) x %>% filter(id %in% y$id),
            epiMatr, intervals)

#10) save matrix in two files:one with only the class profiles and one with all the information from intervals.txt added
dir.create(outDir)
outFiles= split(names(epiMatr), f = names(epiMatr))
Map(function(x,y) write_tsv(x, file = paste(outDir,"/",y,"_epiAnalysis_filtered.txt", sep=""), col_names = T, quote = "none"), 
    epiMatr, outFiles)
