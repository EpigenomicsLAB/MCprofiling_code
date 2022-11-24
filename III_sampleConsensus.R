# This script imports MC profiles from multiple samples and compute a consensus profile
# It first selects common epiloci among all samples
# Then, it compare the profiles between all sample pairs
# It retains only epiloci above the threshold in all sample pairs
# Finally, it computes the consensus MC profile by averaging the relative abundance of each MC among the samples
####################################################################################
library(tidyverse)
library(combinat)
library(parallel)
####Parameters
wd="."
setwd(wd)
inputDir="out_II"
# set number of cores 
n_cores=10
# input files
inputFiles= list.files(inputDir) 
sample_names= sapply(strsplit(inputFiles, split="_"), "[[",1)
inputFiles=split(inputFiles, f=sample_names)
# set output directory
outDir=paste0(wd, "/out_III")
##### Custom functions
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log2(x/y), na.rm = T)

compare_regions=function(df1,df2,pair)
{
  df1=split(df1[-1], f = df1$id)
  df2=split(df2[-1], f = df2$id)
  jsd=tibble("id"= names(df1))
  jsd=Map(function(x,y) JSD(x, y),
              df1, df2)
  jsd=tibble("id"= names(jsd),
             "jsd"=unlist(jsd))
  names(jsd)[2]=pair
  
  return(jsd)
}

compute_avg= function(region,sample.profiles)
{
  reg_classes=Reduce(bind_rows,map(sample.profiles, function(x) x %>% filter(id==region) %>% select(-id)))
  reg_avg=colMeans(reg_classes)
  
  return(reg_avg)
}
#### main
# 1) read input files
classProfiles=map(inputFiles, function(x) read_tsv(paste(inputDir,x,sep="/")))
# 2) select common regions
common_regions= Reduce(intersect, map(classProfiles, function(x) x %>% select(id) %>% pull() %>% return()))
classProfiles = map(classProfiles, 
                    function(x) x %>% filter(id %in% common_regions) %>% arrange(id) %>% return())
# 3) compute JSD distance among all sample pairs
comb=combn(names(classProfiles),2, simplify = F)
jsd=map(comb, function(x) compare_regions(classProfiles[[x[1]]], classProfiles[[x[2]]], paste(x, collapse = "_")))
jsd=Reduce(function(...) inner_join(..., by="id"), jsd)
# 4) save jsd
if(!dir.exists(outDir))
{
  dir.create(outDir)
}
write_tsv(jsd, paste0(outDir,"/jsd_pairs.tsv"), col_names = T, quote = "none")
# 5) select regions with stable MC profiles among samples
jsd = jsd %>% 
  filter(if_any(starts_with("sample"), function(x) x<=0.26))
classProfiles = map(classProfiles, function(x) x %>% filter(id %in% jsd$id) %>% arrange(id) %>% return())
# 6) compute average Profiles
avg_profiles=mclapply(jsd$id, function(x) compute_avg(x, classProfiles), mc.cores = n_cores)
avg_profiles=Reduce(bind_rows, avg_profiles)
jsd=bind_cols(jsd, avg_profiles)
# 7) save consensus profile
write_tsv(jsd, paste(outDir,"/consensus_MCprofiles.txt", sep=""), col_names = T, quote = "none")
