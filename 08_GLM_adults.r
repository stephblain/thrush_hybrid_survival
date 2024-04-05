#############################################################
#Make some plots to look at ancestry and heterozygosity distributions
#Then fit adult survival model
#############################################################



setwd("C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/")


#load packages
#pkgs<-c("tidyverse","HIest","viridis","ggpubr","fields","gsg","mgcv","survival","survminer","mgcViz")
pkgs<-c("tidyverse","ggplot2","HIest","viridis","ggpubr","fields","mgcv","survival","survminer","mgcViz")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

#read in survival functions
source("../thrush_survival/00_functions_survival.R")

#read in file with HIest estimates, individual names, and survival data
HI_thrush<-read.csv("thrush_survival_240106.csv"); HI_thrush$X<-NULL
mean(HI_thrush$ancestry,na.rm=T); 1-mean(HI_thrush$ancestry,na.rm=T)

poly.all<-data.frame(X.co=rep(c(0,0.5,1)), #make triangle polygon for plotting
                     Y.co=rep(c(0,1,0)))
					 
#############################################################
## Check distributions of ancestry across space
#############################################################

#Summary stats by release site - adults and juvies combined

#recode sites by their ancestry type
HI_thrush<-HI_thrush%>%
  filter(!is.na(release_site))%>%
  filter(!is.na(ancestry))%>%
  mutate(release_site=str_trim(release_site))%>%
  mutate(release_site=str_replace(release_site," ","_"))%>%
  mutate(pop_type=dplyr::recode(release_site,"Porpoise_Bay"="coastal","Pacific_Spirit"="coastal",
                                "Kamloops"="inland","Kelowna"="inland",
                                "Williams_Lake"="inland","Tatlayoko"="inland",
                                "Alaska"="hybrid","Bella_Coola"="hybrid","Hope"="hybrid",
                                "Pemberton"="hybrid","Washington"="hybrid"))
								
#look at summary stats by release site

ggplot(data=HI_thrush%>%filter(release_site=="Pemberton"),
       aes(x=as.factor(release_year),fill=tag_type))+geom_bar()+
  scale_fill_manual(values=c("aquamarine4","aquamarine3"))+ggtitle("Pemberton")

HI_thrush%>%
  group_by(release_site,pop_type)%>%
  summarise(mean_ancestry=mean(ancestry),mean_heterozygosity=mean(heterozygosity),count=n())

ggplot(data=HI_thrush,aes(x=fct_reorder(release_site,ancestry,median),y=ancestry,colour=tag_type))+
  geom_hline(yintercept=0.5,lty=2)+geom_boxplot()+
  geom_point(aes(group=tag_type),position = position_jitterdodge(jitter.width=0.2))+
  scale_colour_manual(values=c("aquamarine4","aquamarine3"))+
  xlab("release site")


#Compare ancestry and heterozygosity distributions in and out of the hybrid zone

ggplot()+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="black")+
  geom_point(data=HI_thrush%>%filter(pop_type=="hybrid"),aes(x=ancestry,y=heterozygosity),
             size=1,shape=19,stroke=1)+
  geom_point(data=HI_thrush%>%filter(pop_type%in%c("coastal","inland")),
             aes(x=ancestry,y=heterozygosity,
                 colour=fct_reorder(release_site,ancestry)),size=4,shape=19,stroke=1)+
  coord_equal()+scale_colour_viridis(discrete=T,name="Population",option="C",begin=0.2)+
  xlab("Ancestry")+ylab("Heterozygosity")+
  theme_classic()

#Histograms of ancestry across sites and across longitudes in Pemberton

A<-ggplot()+
  geom_histogram(data=HI_thrush%>%filter(!is.na(release_site)),
                 aes(x=ancestry,fill=fct_reorder(release_site,ancestry)))+
  scale_fill_viridis(discrete=T,option="C",name="release site")+
  facet_wrap(vars(pop_type))+ggtitle("All sites")

B<-ggplot()+
  geom_histogram(data=HI_thrush%>%filter(pop_type%in%c("hybrid")&release_site=="Pemberton"&
                                           !is.na(release_gps.w)),
                 aes(x=ancestry,fill=as.factor(round(release_gps.w,2))))+
  scale_fill_viridis(discrete=T,option="B",name="release GPS W")+
  facet_wrap(vars(tag_type))+ggtitle("Pemberton")

ggarrange(A,B,nrow=2)


## Visualize survival across triangle plots
#Survival to breeding grounds

HI_adult<-HI_thrush%>%filter(age_release%in%c("AHY","SY","ASY")&pop_type=="hybrid")%>%
  mutate(retrieved_archival=as.numeric(retrieved_archival))%>%
  filter(!is.na(retrieved_archival))%>%filter(!is.na(ancestry))
HI_juvie<-HI_thrush%>%filter(age_release=="HY"&pop_type=="hybrid"&release_site=="Pemberton")%>%
  filter(!is.na(ancestry))

p1.tri<-ggplot()+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey50",size=2)+
  geom_point(data=HI_adult,
             aes(x=ancestry,y=heterozygosity,
                 shape=as.factor(retrieved_archival)),colour="grey10",size=2,stroke=1)+
  scale_shape_manual(values=c(1,19),guide="none")+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()


p2.tri<-ggplot()+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey50",size=2)+
  geom_point(data=HI_juvie,
             aes(x=ancestry,y=heterozygosity,
                 shape=as.factor(t5_spring40)),colour="grey10",size=2,stroke=1)+
  scale_shape_manual(values=c(1,19),guide="none")+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()

ggarrange(p1.tri,p2.tri,labels=c("adult","juvenile"))

#Make table for supp mat with all birds in the study

HI_summary<-rbind(HI_juvie,HI_adult)%>%
  mutate(reference=substr(name_in_vcf,1,7),
         latitude=round(release_gps.n,3),longitude=round(release_gps.w,3))%>%
  dplyr::select(tag_type,reference,release_site,latitude,longitude,release_year)%>%
  rename(year=release_year)
  
#write.csv(HI_summary,"../figures/individuals_table.csv")


#############################################################
#Run GLMs to estimate variation in adult survival with ancestry & heterozygosity
#############################################################

HI_adult$release_year_f<-as.factor(paste("y",HI_adult$release_year,sep=""))
HI_adult$sex_binary_f<-as.factor(paste("s",HI_adult$sex_binary,sep=""))

sum(HI_adult$retrieved_archival)

library(car)
library(lme4)
ma1<-glmer(retrieved_archival~ancestry*heterozygosity+
           (1|release_year_f)+(1|sex_binary_f),family=binomial(link="logit"),
    data=HI_adult)
Anova(ma1)
summary(ma1)

HI_adult$fitted_survival<-plogis(predict(ma1))

HI_adult%>%filter(release_year_f=="y2022")%>%pull(retrieved_archival)
#write.csv(HI_adult,"thrush.adult.fitted.survival.csv",row.names=F)
