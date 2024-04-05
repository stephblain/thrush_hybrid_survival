library(tidyverse)
library(Matrix)
library(marked)

thrush1<-read.csv("thrush_survival_ch_240104.csv")

today<-gsub("-","",Sys.Date())

thrush<-thrush1%>%filter(tag_type=="radio")%>%
  mutate(ch=gsub("days_","",ch))%>%
  mutate(year_f=as.factor(paste("f",release_year,sep="")),
         sex_f=as.factor(paste("f",sex_binary,sep="")))%>%
  #filter(!year_f%in%c("f2019","f2020"))%>%
  select(ch,ancestry,heterozygosity,year_f,sex_f)

model.parameters=list(Phi=list(formula=~ancestry*heterozygosity),
                      p=list(formula=~time+year_f+sex_f+ancestry+heterozygosity))

#thrush$ch<-substr(thrush$ch,1,10)

model1="hmmCJS"

cjsfit=crm(thrush,model=model1,hessian=T,
            model.parameters=model.parameters)

sink(paste(model1,"beta",today,"txt",sep="."))
cjsfit$results
sink(file = NULL)

fit1<-predict(cjsfit)
phi1<-fit1$Phi
phi1$est.adj<-phi1$estimate^299
p1<-fit1$p

write.csv(phi1,paste(model1,"fit.phi",today,"csv",sep="."),row.names=F)
write.csv(p1,paste(model1,"fit.p",today,"csv",sep="."),row.names=F)