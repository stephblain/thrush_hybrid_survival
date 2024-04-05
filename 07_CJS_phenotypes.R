library(tidyverse)
library(Matrix)
library(marked)

thrush1<-read.csv("thrush_phenotypes_ch_231130.csv")

today<-gsub("-","",Sys.Date())


########disruptive

thrush<-thrush1%>%filter(tag_type=="radio")%>%
  mutate(ch=gsub("days_","",ch))%>%
  mutate(year_f=as.factor(paste("f",release_year,sep="")),
         sex_f=as.factor(paste("f",sex_binary,sep="")),
         WingShapeTailQuad=WingShapeTail^2,
         BearingTarsusQuad=BearingTarsus^2,
         DepartureDayQuad=DepartureDay^2)%>%
  select(ch,ancestry,heterozygosity,year_f,sex_f,BearingTarsus,
         DepartureDay,WingShapeTail,BearingTarsusQuad,
         DepartureDayQuad,WingShapeTailQuad)

# 
# thrush<-thrush[sample(1:255,30),]%>%
#   mutate(ch=substr(ch,1,20))

model1="hmmCJS"
variable1="BearingTarsus"
model.parameters=list(Phi=list(formula=~BearingTarsus+BearingTarsusQuad),
                      p=list(formula=~time+year_f+sex_f+ancestry+heterozygosity))



cjsfit=crm(thrush,model=model1,hessian=T,
           model.parameters=model.parameters)

sink(paste(model1,variable1,"beta",today,"txt",sep="."))
cjsfit$results
sink(file = NULL)

fit1<-predict(cjsfit)
phi1<-fit1$Phi
phi1$est.adj<-phi1$estimate^299
p1<-fit1$p

write.csv(phi1,paste(model1,variable1,"fit.phi",today,"csv",sep="."),row.names=F)
write.csv(p1,paste(model1,variable1,"fit.p",today,"csv",sep="."),row.names=F)

########mismatch

thrush<-thrush1%>%filter(tag_type=="radio")%>%
  mutate(ch=gsub("days_","",ch))%>%
  mutate(year_f=as.factor(paste("f",release_year,sep="")),
         sex_f=as.factor(paste("f",sex_binary,sep="")))%>%
  rename(mismatch_all=mismatch.all)%>%
  #filter(!year_f%in%c("f2019","f2020"))%>%
  select(ch,ancestry,heterozygosity,year_f,sex_f,BearingTarsus_DepartureDay,
         WingShapeTail_BearingTarsus,WingShapeTail_DepartureDay,mismatch_all)

model1="hmmCJS"
variable1="mismatch_all"
model.parameters=list(Phi=list(formula=~mismatch_all),
                      p=list(formula=~time+year_f+sex_f+ancestry+heterozygosity))



cjsfit=crm(thrush,model=model1,hessian=T,
           model.parameters=model.parameters)

sink(paste(model1,variable1,"beta",today,"txt",sep="."))
cjsfit$results
sink(file = NULL)

fit1<-predict(cjsfit)
phi1<-fit1$Phi
phi1$est.adj<-phi1$estimate^299
p1<-fit1$p

write.csv(phi1,paste(model1,variable1,"fit.phi",today,"csv",sep="."),row.names=F)
write.csv(p1,paste(model1,variable1,"fit.p",today,"csv",sep="."),row.names=F)