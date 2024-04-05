#load packages
pkgs<-c("tidyverse","HIest","viridis","ggpubr","fields")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

setwd("C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/")

#read in survival functions
source("../thrush_survival/00_functions_survival.R")

#read in file with HIest estimates, individual names, and survival data
HI_thrush<-read.csv("thrush_survival_240106.csv"); HI_thrush$X<-NULL
HI_adult<-read.csv("thrush.adult.fitted.survival.csv")
phi1<-read.csv("outputs/hmmCJS.fit.phi.20240104.csv")
p1<-read.csv("outputs/hmmCJS.fit.p.20240104.csv")


poly.all<-data.frame(X.co=rep(c(0,0.5,1)), #make triangle polygon for plotting
                     Y.co=rep(c(0,1,0)))

#recode sites by their ancestry type
HI_thrush<-HI_thrush%>%
  filter(!is.na(release_site))%>%
  mutate(release_site=str_trim(release_site))%>%
  mutate(release_site=str_replace(release_site," ","_"))%>%
  mutate(pop_type=dplyr::recode(release_site,"Porpoise_Bay"="coastal","Pacific_Spirit"="coastal",
                                "Kamloops"="inland","Kelowna"="inland",
                                "Williams_Lake"="inland","Tatlayoko"="inland",
                                "Alaska"="hybrid","Bella_Coola"="hybrid","Hope"="hybrid",
                                "Pemberton"="hybrid","Washington"="hybrid"))

HI_juvie<-HI_thrush%>%filter(age_release=="HY"&pop_type=="hybrid"&release_site=="Pemberton")


##################################################################
#Make Figure 2
##################################################################


cb=0
v.opt="mako"
vir.min.adult<-(min(HI_adult$fitted_survival)-min(phi1$est.adj))/(max(HI_adult$fitted_survival)-min(phi1$est.adj))*(1-cb)+cb
vir.max.juvie<-(max(phi1$est.adj)-min(phi1$est.adj))/(max(HI_adult$fitted_survival)-min(phi1$est.adj))*(1-cb)+cb

min(HI_adult$fitted_survival)
max(phi1$est.adj)

p2.tri<-ggplot(phi1,aes(x=ancestry,y=heterozygosity,colour=est.adj))+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey30",size=0.8)+
  geom_jitter(size=2.5)+
  scale_colour_viridis(option=v.opt,begin=cb,end=vir.max.juvie,name="Survival")+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()+
  theme(panel.background=element_rect(colour="black",size=1),
        legend.position="none")

p1.tri<-ggplot(HI_adult,aes(x=ancestry,y=heterozygosity,colour=fitted_survival))+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey30",size=0.8)+
  geom_jitter(size=2.5)+
  scale_colour_viridis(option=v.opt,begin=vir.min.adult,end=1,name="Survival")+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()+
  theme(panel.background=element_rect(colour="black",size=1),
        legend.position="none")

legend.data=data.frame(Zval=c(HI_adult$fitted_survival,phi1$est.adj),
                       Xval=1,Yval=2)
p2.tile<-ggplot()+
  geom_tile(data=legend.data,aes(Xval,Yval,fill=Zval))+
  scale_fill_viridis(name="Survival",option=v.opt,begin=cb,end=1)


p12<-gridExtra::grid.arrange(p1.tri,p2.tri,get_legend(p2.tile),nrow=1,ncol=3,widths=c(2,2,1))  

str(p1)
p1.time<-p1%>%
  mutate(timef=sprintf("%03d",as.numeric(p1$time)))%>%
  #mutate(ancestry_bin=case_when(ancestry>=0.66~"high",ancestry<0.66&ancestry>0.33~"med",
  #                              ancestry<=0.33~"low"))%>%
  group_by(timef)%>%
  summarise(mean_survival=mean(estimate),time=mean(time))%>%
  filter(time!=300)

p3<-ggplot(p1.time,aes(x=time,y=mean_survival))+
  geom_point(colour="grey30",size=1)+
  ylab("Detection\nprobability")+
  xlab("Days since tagging")+
  theme(panel.background=element_rect(colour="black",size=1))

fig2<-ggarrange(p12,p3,nrow=2,heights=c(0.6,0.4),labels=LETTERS)
fig2

# ggsave(filename="../figures/Fig2.png",device="png",
#        plot=fig2,width=15,height=10,units="cm",bg="white")
# 
# ggsave(filename="../figures/Fig2.pdf",device="pdf",
#        plot=fig2,width=15,height=10,units="cm",bg="white")

dev.off()

##################################################################
#Get selection coefficient estimates
##################################################################

ch<-read.csv("thrush_survival_ch_240104.csv")
nrow(ch%>%filter(age_release=="HY"&release_site=="Pemberton"))


coastal<-phi1%>%filter(heterozygosity<0.25&ancestry<0.25)%>%
  summarise(count=n(),mean_survival=mean(est.adj))%>%pull(mean_survival)
inland<-phi1%>%filter(heterozygosity<0.25&ancestry>0.75)%>%
  summarise(count=n(),mean_survival=mean(est.adj))%>%pull(mean_survival)

high_H<-phi1%>%filter(heterozygosity>0.75)%>%
  summarise(count=n(),mean_survival=mean(est.adj))%>%pull(mean_survival)

F1.s<-1-high_H/mean(coastal,inland)
F1.s

CB_H<-phi1%>%filter(heterozygosity>0.25&heterozygosity<0.75&ancestry<0.4)%>%
  summarise(count=n(),mean_survival=mean(est.adj))%>%pull(mean_survival)
IB_H<-phi1%>%filter(heterozygosity>0.25&heterozygosity<0.75&ancestry>0.6)%>%
  summarise(count=n(),mean_survival=mean(est.adj))%>%pull(mean_survival)

cb.s<-1-CB_H/mean(coastal,inland)
cb.s
ib.s<-1-IB_H/mean(coastal,inland)
ib.s

class_table<-phi1%>%mutate(genomic_class=case_when(heterozygosity<0.25&ancestry<0.25~"coastal",
                                      heterozygosity<0.25&ancestry>0.75~"inland",
                                      heterozygosity>0.75~"F1",
                                      heterozygosity>0.25&heterozygosity<0.75&ancestry<0.4~"coastal cross",
                                      heterozygosity>0.25&heterozygosity<0.75&ancestry>0.6~"inland cross",
                                      TRUE~"other"))%>%
  group_by(genomic_class)%>%
  summarise(count=n(),mean_survival=round(mean(est.adj),3),mean_S=mean(ancestry))%>%
  mutate(age="juvenile")


class_table_2<-phi1%>%mutate(genomic_class=case_when(heterozygosity<0.25&ancestry<0.25~"coastal",
                                      heterozygosity<0.25&ancestry>0.75~"inland",
                                      heterozygosity>0.7~"F1",
                                      heterozygosity>0.25&heterozygosity<0.7&ancestry<0.4~"coastal cross",
                                      heterozygosity>0.25&heterozygosity<0.7&ancestry>0.6~"inland cross",
                                      TRUE~"other"))%>%
  group_by(genomic_class)%>%
  summarise(count=n(),mean_survival=round(mean(est.adj),3),mean_S=mean(ancestry))%>%
  mutate(age="juvenile")

#Selection coefficient estimates with adjusted genomic classes - for supp mat

1-class_table_2%>%filter(genomic_class=="F1")%>%pull(mean_survival)/
  mean(class_table_2%>%filter(genomic_class=="coastal")%>%pull(mean_survival),
       class_table_2%>%filter(genomic_class=="inland")%>%pull(mean_survival))
1-class_table_2%>%filter(genomic_class=="coastal cross")%>%pull(mean_survival)/
  mean(class_table_2%>%filter(genomic_class=="coastal")%>%pull(mean_survival),
       class_table_2%>%filter(genomic_class=="inland")%>%pull(mean_survival))
1-class_table_2%>%filter(genomic_class=="inland cross")%>%pull(mean_survival)/
  mean(class_table_2%>%filter(genomic_class=="coastal")%>%pull(mean_survival),
       class_table_2%>%filter(genomic_class=="inland")%>%pull(mean_survival))


class_table_3<-phi1%>%mutate(genomic_class=case_when(heterozygosity<0.25&ancestry<0.25~"coastal",
                                                     heterozygosity<0.25&ancestry>0.75~"inland",
                                                     heterozygosity>0.8~"F1",
                                                     heterozygosity>0.25&heterozygosity<0.8&ancestry<0.4~"coastal cross",
                                                     heterozygosity>0.25&heterozygosity<0.8&ancestry>0.6~"inland cross",
                                                     TRUE~"other"))%>%
  group_by(genomic_class)%>%
  summarise(count=n(),mean_survival=round(mean(est.adj),3),mean_S=mean(ancestry))%>%
  mutate(age="juvenile")

#Selection coefficient estimates with adjusted genomic classes - for supp mat

1-class_table_3%>%filter(genomic_class=="F1")%>%pull(mean_survival)/
  mean(class_table_3%>%filter(genomic_class=="coastal")%>%pull(mean_survival),
       class_table_3%>%filter(genomic_class=="inland")%>%pull(mean_survival))
1-class_table_3%>%filter(genomic_class=="coastal cross")%>%pull(mean_survival)/
  mean(class_table_3%>%filter(genomic_class=="coastal")%>%pull(mean_survival),
       class_table_3%>%filter(genomic_class=="inland")%>%pull(mean_survival))
1-class_table_3%>%filter(genomic_class=="inland cross")%>%pull(mean_survival)/
  mean(class_table_3%>%filter(genomic_class=="coastal")%>%pull(mean_survival),
       class_table_3%>%filter(genomic_class=="inland")%>%pull(mean_survival))




HI_adult<-HI_thrush%>%filter(age_release%in%c("AHY","SY","ASY")&pop_type=="hybrid")%>%
  mutate(retrieved_archival=as.numeric(retrieved_archival))%>%
  filter(!is.na(retrieved_archival))%>%filter(!is.na(ancestry))
class_table_adults<-HI_adult%>%mutate(genomic_class=case_when(heterozygosity<0.25&ancestry<0.25~"coastal",
                                      heterozygosity<0.25&ancestry>0.75~"inland",
                                      heterozygosity>0.75~"F1",
                                      heterozygosity>0.25&heterozygosity<0.75&ancestry<0.4~"coastal cross",
                                      heterozygosity>0.25&heterozygosity<0.75&ancestry>0.6~"inland cross",
                                      TRUE~"other"))%>%
  group_by(genomic_class)%>%
  summarise(count=n(),sum_survival=round(sum(retrieved_archival),3),mean_S=mean(ancestry))%>%
  mutate(age="adult")%>%
  arrange(age,mean_S)%>%
  select(genomic_class,count,sum_survival)

class_table_adults<-rbind(class_table,class_table_adults)%>%
  arrange(age,mean_S)%>%
  select(age,genomic_class,count,mean_survival)

write.csv(class_table,"../figures/tableS2.csv",row.names=F)

write.csv(class_table_adults,"../figures/tableS3.csv",row.names=F)
