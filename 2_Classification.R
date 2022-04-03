#############################################################################################################
#############################################################################################################
# This script performs a supervised classification of epiloci according to the MC profile                   #
# For each epilocus, the script computes the JSD distance in respect to 5 reference distributions           #
# An epilocus is labelled according to the more similar reference distribution                              #
# In this way, all epiloci are actually classified into 5 groups                                            #
#############################################################################################################
#############################################################################################################

#######input variables
wd = "/home/DatiLab2/antonella/RRBS/DNMTs/out_new//epitabs/"
### Reference distributions
distr = list("distr1"=c(0.4,0.1,0,0.1,0.4),
             "distr2" = c(0.8,0.2,0,0,0),
             "distr3" = c(0,0,0,0.2,0.8),
             "distr4" = c(0,0.25,0.5,0.25,0),
             "distr5" = c(0.2,0.2,0.2,0.2,0.2))
## file containing profiles to be classified
allInfo_file= "WT_3" 
## columns index of MC freq in allinfo file
MC_index = c(15:19)
#######custom functions
KLD <- function(x,y) sum(x * log2(x/y), na.rm = T) #compute Kullback-Libler distance between two distributions
JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2)) #compute Jensen-Shannon distance between two distributions
JSD_archetype= function(df, archetype) apply(df, 1, function(x) JSD(as.vector(t(x)), archetype))

classify_regions=function(df, distr) 
{
  regions=rownames(df)
  distances=lapply(distr, function(x) JSD_archetype(df,x))
  distances=as.data.frame(distances)
  distances$Class = apply(distances, 1, function(x) which.min(x))
  return(distances)
}
#######main
setwd(wd)
allinfo <- read.delim("/home/DatiLab2/antonella/RRBS/DNMTs/out/epitabs/DKO_2_allInfo.txt", stringsAsFactors=FALSE)
# allinfo = read.csv(allInfo_file, sep="", stringsAsFactors=FALSE)
allinfo = WT_3_allInfo
distr_df = allinfo[,MC_index]
rownames(distr_df)=allinfo$id

distances=classify_regions(distr_df, distr) #classify regions
allinfo=merge(allinfo, distances, by.x=1,by.y=0)

### Rimuovi cromosomi sessuali
allinfo = allinfo[-grep("X", allinfo$seqnames),]
allinfo = allinfo[-grep("Y", allinfo$seqnames),]

write.table(allinfo, paste0(gsub(allInfo_file, pattern = "\\.txt", replacement = ""), "_classProfiles.txt"), 
            row.names = F, col.names = T, sep="\t", quote = F) # Save data with classification

################################################################################

allinfo = merge(WT_2_classProfiles, DKO_2_classProfiles, by = "id")

table(allinfo$Class.x, allinfo$Class.y)

### percentuali per riga 
prop.table(table(allinfo$Class.x, allinfo$Class.y), 1)





