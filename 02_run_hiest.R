#module load GCC/11.2.0  OpenMPI/4.1.1 R/4.2.0

#R

source("../HIest.R")
library(tidyverse)


prefix="stitch.ldr01.hwe.s1.fst94.noZ"



G_ld<-read.table(paste(prefix,"hiest_format.012",sep="."),header=F); G_ld$V1<-NULL
inds_ld<-read.table(paste(prefix,"hiest_format.012.indv",sep="."))
loci_ld<-read.table(paste(prefix,"hiest_format.012.pos",sep="."))


rownames(G_ld)<-inds_ld$V1

coastal_s1<-read.table("../coastal.sample1")
inland_s1<-read.table("../inland.sample1")

coastal_s2<-read.table("../coastal.sample2")
inland_s2<-read.table("../inland.sample2")


#replace -1 (missing data) with NA
G_ld[G_ld==-1]<-NA


#extract sample 2 individuals
G_coastal_ld_s2<-G_ld[rownames(G_ld)%in%coastal_s2$V1,]
#colSums(G_coastal_ld_s2,na.rm=T)/(2*colSums(!is.na(G_coastal_ld_s2))) #frequency of alt allele
G_inland_ld_s2<-G_ld[rownames(G_ld)%in%inland_s2$V1,]
#colSums(G_inland_ld_s2,na.rm=T)/(2*colSums(!is.na(G_inland_ld_s2)))

#make table with frequency of alt allele
P_thrush_ld_s2<-data.frame(scaffold=loci_ld$V1,
Locus=1:dim(G_ld)[2],Allele=1,
				P1=(colSums(G_inland_ld_s2,na.rm=T)/(2*colSums(!is.na(G_inland_ld_s2)))),
P2=(colSums(G_coastal_ld_s2,na.rm=T)/(2*colSums(!is.na(G_coastal_ld_s2)))))

#adjust expected frequencies to be more realistic
P_thrush_ld_s2$P1[P_thrush_ld_s2$P1>0.98]<-0.98
P_thrush_ld_s2$P1[P_thrush_ld_s2$P1<0.02]<-0.02
P_thrush_ld_s2$P2[P_thrush_ld_s2$P2>0.98]<-0.98
P_thrush_ld_s2$P2[P_thrush_ld_s2$P2<0.02]<-0.02

write.csv(P_thrush_ld_s2,paste("p_hiest",prefix,"adj98.csv",sep="."),row.names=F)


#########estimate ancestry and heterozygosity


#require 25% of loci for each indidual
inds_ld<-inds_ld[rowSums(!is.na(G_ld))>nrow(loci_ld)*0.25,] #keep names that pass filter
G_ld<-G_ld[rowSums(!is.na(G_ld))>nrow(loci_ld)*0.25,]
length(inds_ld)==nrow(G_ld)&ncol(G_ld)==nrow(P_thrush_ld_s2) #check

#run HIest
HI_ld<-HIest(G_ld,P_thrush_ld_s2,type="allele.count", method = "SANN", iterations = 1000, Cscale = NULL,
   start = c(.5,.5), control = list(fnscale = -1, maxit = 1000))

HI_ld$ind<-inds_ld

write.csv(HI_ld,paste("hiest",prefix,"adj98.csv",sep="."),row.names=F)