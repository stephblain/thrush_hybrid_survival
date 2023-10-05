#for running from cluster: module load R/4.2.0-foss-2021b

#simulate "eggs", adding in some coastals and inlands to account for migration into the hybrid zone


source("/scratch/user/sblain/HIest.R")
library(tidyverse)

setwd("/scratch/user/sblain/thrush_survival/hiest_to2022")

G_thrush<-read.table("stitch_feb2023.filtered25.10kb_noinv.s1.fst95.noz.hiest_format.012",header=F); G_thrush$V1<-NULL
#replace missing data with NA
G_thrush[G_thrush==-1]<-NA

thrush_inds<-read.table("stitch_feb2023.filtered25.10kb_noinv.s1.fst95.noz.hiest_format.012.indv")
P_thrush<-read.csv("P_thrush_feb2023.filtered25.10kb_noinv.s2.fst95.adj98.csv"); P_thrush$X<-NULL



meta_thrush<-read.csv("HI_thrush_230315.csv")
cols_list<-c("name_in_vcf","release_site","release_year","age_release","tag_type","t5_spring40","S")
meta_thrush<-meta_thrush%>%
  mutate(age_release=dplyr::recode(age_release,"HY"="juvie","ASY"="adult","SY"="adult"))%>%
  #remove tag types that don't match age
  filter(!(age_release=="juvie"&tag_type=="archival"))%>%
  filter(!(age_release=="adult"&tag_type=="radio"))%>%
  select(all_of(cols_list))

for(j in 1:10){

meta_adult<-meta_thrush%>%filter(age_release%in%c("juvie")&tag_type=="radio"&
                                   release_site=="Pemberton"&
                                   t5_spring40==1&!is.na(S))
								   
meta_coastal<-meta_thrush%>%filter(release_site=="Pacific Spirit"&
                                   S<0.05&!is.na(S))
meta_coastal<-meta_coastal[sample(1:nrow(meta_coastal),5),] #get 5 coastals
								   
meta_inland<-meta_thrush%>%filter(release_site%in%c("Kamloops","Kamloops ","Kelowna")&
								  S>0.95&!is.na(S))
meta_inland<-meta_inland[sample(1:nrow(meta_inland),10),] #get 10 inlands

meta_adult<-rbind(meta_adult,meta_coastal,meta_inland)

time.id<-gsub("-| |:","",Sys.time()) #get unique time stamp


G_adult<-G_thrush[thrush_inds$V1%in%meta_adult$name_in_vcf,]



G_eggs<-data.frame()


for(i in 1:1000){ #number of juvies to sample
  
	#sample 2 individuals
	s1<-sample(1:nrow(G_adult),2,replace=F)
  
	ind1<-G_adult[s1[1],]
	ind2<-G_adult[s1[2],]
  
	#replace heterozygotes with a randomly sampled 0 or 1
	#if either parent has an NA at that site, offspring has an NA
	ind1[ind1==1&!is.na(ind1)]<-rbinom(length(ind1[ind1==1&!is.na(ind1)]),1,0.5)*2
	ind2[ind2==1&!is.na(ind2)]<-rbinom(length(ind2[ind2==1&!is.na(ind2)]),1,0.5)*2
	H1<-(ind1+ind2)/2
  
	G_eggs<-rbind(G_eggs,H1)}

HI_eggs<-HIest(G_eggs,P_thrush,type="allele.count",method = "SANN",
			   iterations = 1000, Cscale = NULL,start = c(.5,.5),
			   control = list(fnscale = -1, maxit = 1000))

HI_eggs<-HI_eggs%>%select(S,H)%>%mutate("run"=paste("run",time.id,sep=""))




write.csv(HI_eggs,paste("eggs","hiest",time.id,"csv",sep="."))

print(j)
}
