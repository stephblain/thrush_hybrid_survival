##Run multistate CJS model in marked

library(tidyverse)
library(Matrix)
library(marked)

model1="hmmMSCJS"

thrush1<-read.csv("thrush_survival_ch_240104.csv")

today<-gsub("-","",Sys.Date())

prepare_admb=function()
{
  Sys.setenv(PATH = paste("/scratch/user/sblain/tools/admb/bin;/scratch/user/sblain/tools/admb/utilities;",
                          Sys.getenv("PATH"), sep = ";"))
  Sys.setenv(ADMB_HOME = "/scratch/user/sblain/tools/admb")
  invisible()
}

prepare_admb()

thrush.m<-thrush1%>%filter(tag_type=="radio")%>%
  mutate(year_f=as.factor(paste("f",release_year,sep="")),
         sex_f=as.factor(paste("f",sex_binary,sep="")))%>%
  select(-ch)%>%rename(ch=ch_multi)%>%
  select(ch,ancestry,heterozygosity,year_f,sex_f)

model.parameters=list(S=list(formula=~ancestry*heterozygosity),
                      Psi=list(formula=~ancestry*heterozygosity),
                      p=list(formula=~time+year_f+sex_f+ancestry+heterozygosity))

cjsfit=crm(thrush.m,model=model1,hessian=T,
           model.parameters=model.parameters,
           strata.labels=c("A","B","C"))

sink(paste(model1,"beta",today,"txt",sep="."))
cjsfit$results
sink(file = NULL)

fit1<-predict(cjsfit)
S1<-fit1$S
psi1<-fit1$Psi
p1<-fit1$p

write.csv(S1,paste(model1,"fit.s",today,"csv",sep="."),row.names=F)
write.csv(psi1,paste(model1,"fit.psi",today,"csv",sep="."),row.names=F)
write.csv(p1,paste(model1,"fit.p",today,"csv",sep="."),row.names=F)