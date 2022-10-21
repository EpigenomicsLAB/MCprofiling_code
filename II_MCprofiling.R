# This script computes MC profiles from the epiStatProfiler output (epiAnalysis.txt file)
# The scripts can work on several input files 
####################################################################################
library(tidyverse)
library(gtools)
## Main variables
wd="."
setwd(wd)
inputDir="out_I"
# import the epiAnalysis.txtfiles from inputDir
matrixFiles=list.files(inputDir) %>%
 .[matches(vars = ., match = "epiAnalysis")]
sample_names= sapply(strsplit(matrixFiles, split="_"), "[[",1)
# generate a named list of files
matrixFiles=split(matrixFiles, f=sample_names)
# set output dir
outDir=paste0(wd, "/out_II")
####Custom functions
compute_classProfiles=function(sampleEpi)
{
  epialleles=as_tibble(permutations(2,4,c(0,1), repeats.allowed = T))
  epialleles=mutate(epialleles, class=as.character(rowSums(epialleles)))
  epialleles=unite(epialleles, col = "epi", c("V1":"V4"), remove = T, sep = "")
  sampleEpi=left_join(sampleEpi, epialleles, by=c("Var1"="epi"))
  sampleEpi=sampleEpi[c(2,3,5)]
  sampleEpi=sampleEpi %>% 
    dplyr::group_by(id,class) %>% 
    dplyr::summarise_at(vars(Freq),list(Freq = sum))
  sampleEpi=spread(sampleEpi, key = class, value = Freq)
  sampleEpi[is.na(sampleEpi)]=0
  sampleEpi[-1]=sampleEpi[-1]/rowSums(sampleEpi[-1])
  
  return(as_tibble(sampleEpi))
}
#######main
#1) read input files
epiMatr=map(matrixFiles, 
            function(x) read_tsv(paste(inputDir,x,sep="/")))
#2) compute class distributions
classProfiles=map(epiMatr, function(x) suppressWarnings(compute_classProfiles(x)))

#3) save matrix in two files:one with only the class profiles and one with all the information from intervals.txt added
if(!dir.exists(outDir))
{
  dir.create(outDir)
}

outFiles= split(names(classProfiles), f = names(classProfiles))
Map(function(x,y) write_tsv(x, path = paste(outDir,"/",y,"_MCprofiles.txt", sep=""), col_names = T, quote = "none"),
    classProfiles, outFiles)
