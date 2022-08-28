# This script assignes an MC profile to the nearest methylation pattern (MP)
# First, it computes the distance from 5 reference profiles
# Then, it assigns the profile to the MP corresponding to the most similar reference profile 
# The present script works on the consensus MC profile, but can be arranged to work on individual samples' 
# MC profiles (out_II) by changing the value of MC_index
####################################################################################
library(tidyverse)
######
####Parameters
wd="."
setwd(wd)
inputDir="out_III"
# input files
inputFile= "consensus_MCprofiles.txt"
# MC index (columns indicating the relative abundance of MCs)
MC_index=c(5:9)
# set output directory
outDir=paste0(wd, "/out_IV")
##### reference profiles
distr = list("D1"=c(0.8,0.2,0,0,0),
             "D2" = c(0,0,0,0.2,0.8),
             "D3" = c(0.4,0.1,0,0.1,0.4),
             "D4" = c(0,0.25,0.5,0.25,0),
             "D5" = c(0.2,0.2,0.2,0.2,0.2))
##### custom functions 
KLD <- function(x,y) sum(x * log2(x/y), na.rm = T) #compute Kullback-Libler distance between two distributions
JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2)) #compute Jensen-Shannon distance between two distributions
JSD_archetype= function(df, archetype) apply(df, 1, function(x) JSD(as.vector(t(x)), archetype))

classify_regions=function(df, distr) 
{
  regions=rownames(df)
  distances=map(distr, function(x) JSD_archetype(df,x))
  distances=as.data.frame(distances)
  distances$Class = apply(distances, 1, function(x) names(distances)[which.min(x)])
  distances$id=rownames(distances)
  distances=as_tibble(distances) 
  
  return(distances)
}

##### main
# 1) import MC profiles to be classified
profiles=read_tsv(paste(inputDir,inputFile,sep="/"))
# 2) compute the distances from the reference profiles and assign to MPs
distr_df = as.data.frame(profiles[,MC_index])
rownames(distr_df)=profiles$id
distances=classify_regions(distr_df, distr)
allinfo=left_join(profiles, distances, by="id")
# 3) save output
dir.create(outDir)
write_tsv(allinfo, paste(outDir,"/MP_classification.txt", sep=""), col_names = T, quote = "none")


