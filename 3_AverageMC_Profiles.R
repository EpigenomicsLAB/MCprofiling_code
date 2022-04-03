#############################################################################################################
#############################################################################################################
# This script calculates average MC profiles from diffirent replicates using as input the files             #
# coming from EpiStatProfiler.                                                                             #
# New tables containing all the information will be saved in the working directory.                                                                                                          #
#                                                                                                           #
#############################################################################################################
#############################################################################################################
library(rlist)
wd="/home/DatiLab2/Condivisa_2021/Samples_DNMTsKO/"
setwd(wd)
group="D1KO"
files=list.files(wd) %>% .[matches(paste0(group,".*classProfiles.txt"), vars=.)]
profiles=lapply(files[c(1,2)], function(x) read.csv(x, header = F, sep = "\t"))
profiles2=lapply(files[c(3)], function(x) read.csv(x, header = T, sep = "\t"))
profiles=append(profiles,profiles2)
profiles=lapply(profiles,"[",,c(1, 15:19))
#### V1 renamed as id; V2:V6 renames as sample_<class>
names=c("id", paste0("X", seq(0,4)))
profiles=lapply(profiles, function(x) {names(x)=names; return(x)})
regions=lapply(profiles,"[",,"id")
regions=lapply(regions, function(x) {x=as.data.frame(x); x$sample="TRUE"; return(x)})
regions=Reduce(function(...) merge(..., all=T,by=1), regions)
names(regions)=c("id",paste0("S", seq(length(profiles))))
regions[is.na(regions)]="FALSE"
regions$count=apply(regions,1, function(x) sum(x=="TRUE"))

#index= lapply(seq(nrow(regions)), function(x) which(regions[x,-1]=="TRUE"))
#unlist(mclapply(seq(nrow(regions))[c(1:10)], function(x) compute_jsd(as.character(allsamples_regions[x,1]),profiles[index[[x]]]), mc.cores = 10))
#compute average profile

compute_avg= function(region,sample.profiles)
{
  reg_classes=Reduce(rbind,lapply(sample.profiles, function(x) x %>% dplyr::filter(id==region) %>% dplyr::select(-id)))
  reg_avg=as.data.frame(t(colMeans(reg_classes)))
  return(reg_avg)
}


avg_profiles=mclapply(as.character(regions$id), function(x) compute_avg(x, profiles), mc.cores = 10)

regions=cbind(regions,Reduce(rbind, avg_profiles))

#write.table(regions,"/home/DatiLab2/Giulia/ampliconi2.0/rep_variability_092021/regions_all.txt", col.names = T, row.names = F, quote = F, sep="\t")

compute_jsd= function(region, sample.profiles)
{
  reg_classes=Reduce(rbind,lapply(sample.profiles, function(x) x %>% dplyr::filter(id==region) %>% dplyr::select(-id)))
  reg_avg=colMeans(reg_classes)
  KLD <- function(x,y) sum(x * log2(x/y), na.rm = T)
  distances=apply(reg_classes,1, function(x) KLD(x,reg_avg))
  jsd=sqrt(sum(distances*(1/length(distances))))
  
  return(jsd)
}

regions$jsd=unlist(mclapply(as.character(regions[,1]), function(x) compute_jsd(x,profiles), mc.cores = 10))

write.table(regions,paste0("/home/DatiLab2/Giulia/ampliconi2.0/rep_variability_092021/regions_all_",group, ".txt"), col.names = T, row.names = F, quote = F, sep="\t")

######### merge 2 tables

wd="/home/DatiLab2/Giulia/ampliconi2.0/rep_variability_092021/"
setwd(wd)
groups=c("WT","DKO")
##import table 1
tab1=list.files(wd) %>% .[matches(groups[1], vars=.)] %>% .[2]
tab1=read.csv(tab1, header = T, sep = "\t")
names(tab1)[-1]=paste(groups[1], names(tab1)[-1], sep="_")
tab1=tab1[tab1$WT_count>2,]
##import table 2
tab2=list.files(wd) %>% .[matches(groups[2], vars=.)] %>% .[2]
tab2=read.csv(tab2, header = T, sep = "\t")
names(tab2)[-1]=paste(groups[2], names(tab2)[-1], sep="_")
tab2=tab2[tab2$DKO_count>2,]
## merge tables
tab=inner_join(tab1, tab2, by="id")
## calculate the third shannon and add it to the table
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log2(x/y), na.rm = T)
tab$comp_shannon=sapply(seq(nrow(tab)), function(x) JSD(as.vector(tab[x,c(8:12)]), as.vector(tab[x,c(26:30)])))
tab$comp_shannon[is.nan(tab$comp_shannon)]=0
## save output
#write.table(tab,paste0(paste(groups, collapse="_") ,".txt"), col.names = T, row.names = F, quote = F, sep="\t")
tab$WT_se=tab$WT_jsd/tab$WT_count
tab$DKO_se=tab$DKO_jsd/tab$DKO_count
tab$stat=tab$comp_shannon/(sqrt(tab$WT_se+tab$DKO_se))
plot(tab$comp_shannon, tab$stat,cex=0.01)
# 1) plot complessivo con soglie
plot(tab$comp_shannon, tab$DKO_jsd, cex=0.1)
plot(tab$comp_shannon, tab$WT_jsd, cex=0.1)
# 2) taglio la tabella
vars=tab %>% filter(comp_shannon >0.7 & DKO_jsd<=0.4 & WT_jsd<=0.4)
vars=tab %>% filter(stat>1.5)
hist(vars$comp_shannon)
# 3) cross table
table(vars$WT_avg_Class, vars$DKO_avg_Class)
###########test whether old function and new function for jsd work in the same way
# tab <- read.delim("/home/DatiLab2/Giulia/ampliconi2.0/rep_variability_092021/WT_DKO.txt", stringsAsFactors=FALSE)
# length(which(tab$WT_count>2 & tab$DKO_count>2))
# length(which(tab$WT_count>3 & tab$DKO_count>2))
# length(which(tab$WT_count>4 & tab$DKO_count>2))
# tab_all=tab[tab$WT_count>4 & tab$DKO_count>2,]
# 
# JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
# KLD <- function(x,y) sum(x * log2(x/y), na.rm = T)
# 
# a=as.vector(t(tab_all[50,c(8:12)]))
# b=as.vector(t(tab_all[50,c(18:22)]))
# JSD(a,b)
# a=tab_all[50,c(8:12)]
# names(a)=as.character(seq(0,4))
# b=tab_all[50,c(18:22)]
# names(b)=as.character(seq(0,4))
# reg_classes=rbind(a,b)
# reg_avg=colMeans(reg_classes)
# KLD <- function(x,y) sum(x * log2(x/y), na.rm = T)
# distances=apply(reg_classes,1, function(x) KLD(x,reg_avg))
# sqrt(sum(distances*(1/length(distances))))
####classification of average profile and of individual sample profiles
### 1) table of individual profiles
wd="/home/DatiLab2/Condivisa_2021/Samples_DNMTsKO/"
setwd(wd)
group="D1KO"
files=list.files(wd) %>% .[matches(paste0(group,".*classProfiles.txt"), vars=.)]
profiles=lapply(files[c(1,2)], function(x) read.csv(x, header = F, sep = "\t"))
profiles2=lapply(files[c(3)], function(x) read.csv(x, header = T, sep = "\t"))
profiles=append(profiles,profiles2)
lapply(profiles, head)
profiles=lapply(profiles,"[",,c(1,25))
profiles=Reduce(function(...) merge(..., all=T,by=1), profiles)
names(profiles)=c("id",paste(paste0("S", seq(length(profiles)-1)), "Class", sep="_"))
### 2) table of average profiles
wd="/home/DatiLab2/Giulia/ampliconi2.0/rep_variability_092021/"
setwd(wd)
files=list.files(wd) %>% .[matches(paste0("regions_all_",group,".txt"), vars=.)]
tab=read.csv(files, header = T, sep = "\t")
tab=merge(tab, profiles,by="id")
rm(profiles, profiles2)
distr = list("distr1"=c(0.4,0.1,0,0.1,0.4),
             "distr2" = c(0.8,0.2,0,0,0),
             "distr3" = c(0,0,0,0.2,0.8),
             "distr4" = c(0,0.25,0.5,0.25,0),
             "distr5" = c(0.2,0.2,0.2,0.2,0.2))
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
distr_df = tab[c(6:10)]
rownames(distr_df)=tab$id
distances=classify_regions(distr_df, distr) #classify regions
distances=distances[6]
names(distances)="avg_Class"
tab=merge(tab, distances, by.x=1,by.y=0)
tab$conc=apply(tab[c(12:14)],1,function(x) length(unique(x[!is.na(x)])))
tab$conc_avg=apply(tab[c(12:15)],1,function(x) length(which(x[-4][!is.na(x[-4])]==x[4]))/length(x[-4][!is.na(x[-4])]))
write.table(tab,paste0("/home/DatiLab2/Giulia/ampliconi2.0/rep_variability_092021/regions_allInfo_",group, ".txt"), col.names = T, row.names = F, quote = F, sep="\t")
