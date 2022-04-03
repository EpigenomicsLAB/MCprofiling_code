#############################################################################################################
#############################################################################################################
# This script filters out regions from EpiStatProfiler output                                               #
# First, it filters out regions based on their coverage (min-max readNum)                                   #
# It then flags overlapping regions in order to remove them                                                 #
# The scripts can work on several input files                                                               #
#############################################################################################################
#############################################################################################################

####Parameters
wd="/home/DatiLab2/antonella/RRBS/DNMTs/out/"
setwd(wd)

# the input files should be the intervals.txt file created by EpiTool
inputFiles= list("./SRR9016929_sorted_intervals.bed", "./SRR9016930_sorted_intervals.bed",
                 "./SRR9016934_sorted_intervals.bed", "./SRR9016935_sorted_intervals.bed", "./SRR9016943_sorted_intervals.bed")
# the input files should be the epi.txt file created by EpiTool. If more than one file is provided, items should be placed in the same order than inputFiles
matrixFiles= list("./SRR9016929_sorted_epiAnalysis.txt", "./SRR9016930_sorted_epiAnalysis.txt",
                  "./SRR9016934_sorted_epiAnalysis.txt", "./SRR9016935_sorted_epiAnalysis.txt", "./SRR9016943_sorted_epiAnalysis.txt")

#output file names
output=list("WT_1", "WT_2", "D1KO_1", "D1KO_2", "DKO_1")
#minimum coverage required for a region to be retained
min_readNum = 50

#minimum coverage required for a region to be retained
#regions with coverage above the 99.9 percentiles will be fileterd out
#the user can modify the script behavior by setting a value for maxreadNum 
#the custom values have to be provided as a list of the same length of inputFiles
maxreadNum=list()

####Custom functions
removeOverlaps=function(overlap, cores)
{
  rownames(overlap)=seq(nrow(overlap))
  overlap=split(overlap,f = overlap$seqnames)
  overlap=mclapply(overlap, function(x) {
    i=1
    retain=c()
    while(i<=nrow(x))
    {
      retain=append(retain,i)
      chr=x[i,1]
      end=x[i,3]
      df=x[c(i:nrow(x)),]
      Next=min(which(df$seqnames==chr & df$start>end))
      i=i+Next-1
    }
    return(x[retain,])
  },mc.cores = cores)
  return(ldply(overlap, data.frame, .id = c()))
}

compute_classProfiles=function(sampleEpi)
{
  require(gtools)
  epialleles=permutations(2,4,c(0,1), repeats.allowed = T)
  
  epialleles=data.frame("epi"=as.numeric(apply(epialleles,1,function(x) paste(x, collapse = ""))),
                        "class"=as.character(rowSums(epialleles)))
  sampleEpi=merge(sampleEpi, epialleles, by=1)
  sampleEpi=sampleEpi[c(2,3,5)]
  require(dplyr)
  sampleEpi=sampleEpi %>% group_by(class,id) %>% summarise_each(funs(sum))
  require(tidyr)
  sampleEpi=as.data.frame(sampleEpi)
  sampleEpi=spread(sampleEpi, key = id, value = Freq)
  rownames(sampleEpi)=sampleEpi$class
  sampleEpi$class=NULL
  sampleEpi=as.data.frame(t(sampleEpi))
  sampleEpi[is.na(sampleEpi)]=0
  sampleEpi=as.data.frame(apply(sampleEpi,2,function(x) x/rowSums(sampleEpi)))
  return(sampleEpi)
}
#######main
#1) read input files
intervals=lapply(inputFiles, function(x) read.csv(x,sep=""))
#2) assign region id (Chr_start_end)
intervals=lapply(intervals, function(x) {x$id=apply(x[c(1:3)],1,function(x) paste(gsub(x,pattern = " ", replacement = "") ,collapse="_"));return(x)})
#3) filters regions with coverage below minThresh
intervals=lapply(intervals, function(x) return(x[x$num_reads >=50,]))
#4) filters out regions above maxThresh
if(length(maxreadNum)==0)
{
  maxreadNum=lapply(intervals, function(x) quantile(x$num_reads, seq(0,1,0.001))[1000])
}
intervals=lapply(seq(length(intervals)), function(x) return(intervals[[x]][intervals[[x]]$num_reads <= maxreadNum[[x]],]))
#5) remove duplicated rows
intervals=lapply(intervals, function(x) return(x[!duplicated(x$id),]))
lapply(intervals, dim)
#6) select non overlapping regions
intervals2=lapply(intervals, function(x) removeOverlaps(as.data.frame(x), 5))
lapply(intervals2, dim)
#7) load epiMatrix 
epiMatr=lapply(matrixFiles, function(x) read.csv(x, sep=""))
#8) get epialleles of selected regions in intervals
epiMatr=lapply(seq(epiMatr),  function(x) return(epiMatr[[x]][epiMatr[[x]]$id %in% intervals2[[x]]$id,]))
#9) compute class distributions
classProfiles=mclapply(epiMatr, function(x) compute_classProfiles(x), mc.cores=4)
lapply(classProfiles, dim)
#10) save matrix in two files:one with only the class profiles and one with all the information from intervals.txt added
dir.create(paste(wd, 'epitabs', sep='/'))
lapply(seq(length(classProfiles)), function(x) write.table(classProfiles[[x]], file = paste(wd,"/",output[[x]],"_classDist.txt", sep=""), row.names = T, col.names = T, quote = F, sep="\t"))
merged=lapply(seq(length(classProfiles)), function(x) merge(intervals2[[x]], classProfiles[[x]], by.x=14, by.y=0))
lapply(merged,dim)
lapply(seq(length(merged)), function(x) write.table(merged[[x]], file = paste(wd,"/epitabs/",output[[x]],"_allInfo.txt", sep=""), row.names = F, col.names = T, quote = F, sep="\t"))
