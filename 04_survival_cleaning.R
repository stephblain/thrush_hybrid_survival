###################################################################
#Clean up raw radio tag data to estimate survival
###################################################################

#load packages
pkgs<-c("tidyverse","viridis","ggpubr","sf","rnaturalearth")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

setwd("C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/raw/")

#This is the output from Hannah's scripts 1 and 2 - on Delmore github under survival analysis
#see those scripts for additional filtering
lat.df<-read.csv("2.filtered_data_Jul_1.csv")

#read in meta data - output from 01_phenotypes_cleaning.R
meta_thrush<-read.csv("../thrush_meta_240104.csv")%>%
  select(-ancestry,-heterozygosity)

#read in adult recapture data - previously compiled for 2010-2022
# and new list for 2023 recaptures
adult_recap_2023<-read.csv("adults_recaptured_2023.csv")
adult_recap_20102022<-read.csv("survival_data_motus_by_latitude_updated.csv")%>%
  filter(tag_type=="archival"&age_release!="HY")

#load map data
world <- ne_countries(scale = "medium", returnclass = "sf")

hi.noZ<-read.csv("hiest.stitch.ldr01.hwe.s1.fst94.noZ.adj98.csv")%>%
  mutate(S=S*(-1)+1)%>%select(-logLik)%>%
  rename(name_in_vcf=ind,ancestry=S,heterozygosity=H)
hi.noZno4<-read.csv("hiest.stitch.ldr01.hwe.s1.fst94.noZ.no4.adj98.csv")%>%
  mutate(S=S*(-1)+1)%>%select(-logLik)%>%
  rename(name_in_vcf=ind,ancestry_no4=S,heterozygosity_no4=H)
hi.scaf4<-read.csv("hiest.stitch.ldr01.hwe.s1.fst94.noZ.scaf4.adj98.csv")%>%
  mutate(S=S*(-1)+1)%>%select(-logLik)%>%
  rename(name_in_vcf=ind,ancestry_scaf4=S,heterozygosity_scaf4=H)
hi.all<-left_join(left_join(hi.noZ,hi.noZno4),hi.scaf4)


meta_thrush<-left_join(meta_thrush,hi.all)

p_hi<-read.csv("p_hiest.stitch.ldr01.hwe.s1.fst94.noZ.adj98.csv")
nrow(p_hi) #number of sites
#distribution of sites on scaffolds
ggplot(data=p_hi,aes(x=scaffold))+
  geom_bar()+ylab("sites")

poly.all<-data.frame(X.co=rep(c(0,0.5,1)), #make triangle polygon for plotting
                     Y.co=rep(c(0,1,0)))


ggplot()+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="black")+
  geom_point(data=meta_thrush,aes(x=ancestry,y=heterozygosity),size=1,shape=1,stroke=1)+
  coord_equal()+xlab("ancestry")+ylab("heterozygosity")



#clean up detections
#based on visual surveys of routes
lat.df<-lat.df%>%
  
  filter(!is.na(tagDeployStart))%>%
  
  #remove towers with many faulty detections
  filter(!recvDeployName%in%c("LacEdouard-Champs","Lambs Gap")&recvDeployLon<(-70))%>%
  
  #unlikely detections based on time of year and location (north in winter / south in summer)
  filter(!(name_in_vcf=="BI03H07"&(detectdoy>292|detectdoy<291)))%>%
  filter(!name_in_vcf%in%c("AI01H17","AI02H10"))%>%
  filter(!(recvDeployName=="Estero Punta Banda"&detectdoy%in%334:335))%>%
  
  #unlikely detections based on other detections for the same bird
  #ex. inland migration route with one coastal tower - remove coastal tower
  filter(!(recvDeployName=="Golfo de Santa Clara - RV parl"&
             name_in_vcf%in%c("AI01H14","AI02H08","BH31H12","BI03H07","CI03H02_S84_L001","BI02H10")))%>%
  filter(!(name_in_vcf%in%c("AH28H01","AH28H13","AI01H14","AI02H06","AI02H08")&
             recvDeployName=="Drasher"))%>%
  filter(!(name_in_vcf=="CH30H05_S68_L001"&detectdoy>230&detectdoy<234))%>%
  filter(!(name_in_vcf%in%c("BI04H14","BI07H03")&recvDeployName=="Phoenix"))%>%
  filter(!(name_in_vcf%in%c("BH30H01","CH31H03_S77_L001","CH31H05_S79_L001")&
             recvDeployName=="Mackay Island NWR, NC"))%>%
  filter(!(name_in_vcf%in%c("CI03H05_S87_L001")&recvDeployName=="McGill_Bird_Observatory"))%>%
  filter(!(name_in_vcf=="DH31H08"&recvDeployName=="Hopkins Forest"))%>%
  
  #no detections prior to year 2
  filter(!(name_in_vcf%in%c("BI03H05_S181_L001","BI05H03_S185_L001")))%>%
  
  
  #apply general day by location cutoffs
  filter(!(recvDeployLat<30&detectdoy>136&detectdoy<243))%>% #below 30N, June 15 to August 30
  filter(!(recvDeployLat>35&(detectdoy>334|detectdoy<59)))%>% #above 35N, Dec to Feb
  
  mutate(deployDay=substr(tagDeployStart,1,10),tsDay=substr(ts,1,10))%>%
  filter(!is.na(name_in_vcf))%>%
  mutate(tagDays=as.integer(gsub(" days"," ",difftime(tsDay,deployDay,units="days")))) 


lat.sum<-lat.df%>% #next use case when to make column retrievals column and retrieval counts
  mutate(retrieval_time=case_when(recvDeployLat>40&detectdoy>200~"t1_fall40",
                              recvDeployLat>25&recvDeployLat<40&detectdoy>200~"t2_fall25",
                              recvDeployLat<25~"t3_winter00",
                              recvDeployLat>25&recvDeployLat<40&detectdoy<200~"t4_spring25",
                              recvDeployLat>40&detectdoy<200~"t5_spring40"))%>%
  group_by(retrieval_time,name_in_vcf)%>%
  summarise(retrieval_count=n())%>%
  ungroup()%>%mutate(retrieval_count=1)%>%
  pivot_wider(names_from=retrieval_time,values_from=retrieval_count)

lat.sum<-left_join(lat.sum,lat.df%>%group_by(name_in_vcf)%>%
  summarise(taggedDays=max(tagDays)))

###################################################################
#format capture histories for input to CJS models
###################################################################


ch.df<-data.frame()
ch_multi.df<-data.frame()

#for each individual in juvenile survival dataset
#released in Pemberton, not translocated
for(i in meta_thrush%>%filter(tag_type=="radio"&age_release=="HY"&release_site=="Pemberton")%>%
    filter(!name_in_vcf%in%c("BH30H01","BH29H01","BH29H03","BH29H04","BH29H05","BH29H07"))%>%
    pull(name_in_vcf)){
  x1<-rep(0,300)
  x2<-rep(0,300)
  
  y<-lat.df%>%filter(name_in_vcf==i)%>%
    mutate(ch_multi=case_when(recvDeployLat>40&tagDays<150~"A",
                                recvDeployLat<40~"B",
                                recvDeployLat>40&tagDays>150~"C"))%>%
    select(name_in_vcf,tagDays,ch_multi)%>%distinct()
  
  x2[y%>%filter(tagDays<301)%>%pull(tagDays)]<-
    y%>%filter(tagDays<301)%>%pull(ch_multi)
  x2[1]<-"A"
  
  y1<-y%>%pull(tagDays)
  x1[c(1,y1[y1<301])]<-1 #add 1 to detections because bird was caught the day it was tagged
  if(length(y1)>0){if(max(y1)>300){
    x1[300]<-1
    x2[300]<-"C"    }} #fill in last day if bird caught later
  
  ch.df<-rbind(ch.df,c(i,x1))
  ch_multi.df<-rbind(ch_multi.df,c(i,x2))
  
  }

colnames(ch.df)[1]<-"name_in_vcf"
colnames(ch_multi.df)[1]<-"name_in_vcf"

ch.df<-ch.df%>%unite(ch,-name_in_vcf,sep="")
ch_multi.df<-ch_multi.df%>%unite(ch_multi,-name_in_vcf,sep="")
ch.df<-ch.df%>%left_join(ch_multi.df)%>%left_join(meta_thrush)%>%
  mutate(ch=paste("days",ch,sep="_"))%>% #stops capture history from being saved as a very large number
  select(name_in_vcf,ancestry,heterozygosity,
         ancestry_no4,heterozygosity_no4,ancestry_scaf4,heterozygosity_scaf4,
         release_site,release_year,tag_type,sex_binary,age_release,
         release_gps.n,release_gps.w,motustagid,ch,ch_multi)


#write.csv(ch.df,"C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/thrush_survival_ch_240104.csv")

###################################################################
##Get a binary estimate of survival at different time points
#Not used for model-fitting
###################################################################

#fill in previous cols with 1 if detected at a later time point
lat.new<-lat.sum%>%
  mutate(t4_spring25=case_when(t5_spring40==1~1,TRUE~t4_spring25))%>%
  mutate(t3_winter00=case_when(t4_spring25==1~1,TRUE~t3_winter00))%>%
  mutate(t2_fall25=case_when(t3_winter00==1~1,TRUE~t2_fall25))%>%
  mutate(t1_fall40=case_when(t2_fall25==1~1,TRUE~t1_fall40))

#add survival to metadata
HI_thrush<-left_join(meta_thrush,lat.new)%>%
  select(name_in_vcf,ancestry,heterozygosity,
         ancestry_no4,heterozygosity_no4,ancestry_scaf4,heterozygosity_scaf4,
         release_site,release_year,tag_type,sex_binary,age_release,
         release_gps.n,release_gps.w,motustagid,
         t1_fall40,t2_fall25,t3_winter00,t4_spring25,t5_spring40,taggedDays)

#apply filters to only keep birds whose survival we care about
HI_thrush<-HI_thrush%>%
  filter(!tag_type%in%c("captive"))%>% #remove captive birds and birds translocated from whistler
  filter(!name_in_vcf%in%c("BH30H01","BH29H01","BH29H03","BH29H04","BH29H05","BH29H07"))%>%
  #mutate(retrieved_archival=dplyr::recode(tolower(retrieved_archival),"y"=1,"n"=0))%>%
  #remove tag types that don't match age and only keep Pemberton juvies
  filter(!(age_release%in%c("HY")&tag_type%in%c("archival")))%>%
  filter(!(age_release%in%c("SY","ASY")&tag_type%in%c("radio")))%>%
  filter(!(age_release=="HY"&release_site!="Pemberton"))


#fill in zeroes for juvie survival columns
HI_thrush<-HI_thrush%>%
  mutate(taggedDays=as.numeric(taggedDays))%>% #so that not combining num and int later
  mutate(t1_fall40=case_when(age_release=="HY"&is.na(t1_fall40)~0,TRUE~t1_fall40),
         t2_fall25=case_when(age_release=="HY"&is.na(t2_fall25)~0,TRUE~t2_fall25),
         t3_winter00=case_when(age_release=="HY"&is.na(t3_winter00)~0,TRUE~t3_winter00),
         t4_spring25=case_when(age_release=="HY"&is.na(t4_spring25)~0,TRUE~t4_spring25),
         t5_spring40=case_when(age_release=="HY"&is.na(t5_spring40)~0,TRUE~t5_spring40),
         taggedDays=case_when(age_release=="HY"&is.na(taggedDays)~0,TRUE~taggedDays))

#count birds that survived for over a year as survived to 40
HI_thrush<-HI_thrush%>%
  mutate(t1_fall40=case_when(taggedDays>365~1,TRUE~t1_fall40),
         t2_fall25=case_when(taggedDays>365~1,TRUE~t2_fall25),
         t3_winter00=case_when(taggedDays>365~1,TRUE~t3_winter00),
         t4_spring25=case_when(taggedDays>365~1,TRUE~t4_spring25),
         t5_spring40=case_when(taggedDays>365~1,TRUE~t5_spring40))

###################################################################
##Look at individual towers and birds for filtering
#results applied above in lines 62-97
###################################################################


#print data for particular towers
tower.names<-c("Estero Punta Banda","Golfo de Santa Clara - RV parl")
tower.check<-lat.df%>%filter(recvDeployName%in%tower.names)%>%
   mutate(ts=substr(ts,1,10),tagDeployStart=substr(tagDeployStart,1,10))%>%
   group_by(name_in_vcf,tagDeployStart,ts,recvDeployName,recvDeployLat,recvDeployLon,motusTagID)%>%
     summarise(nDetect=n(),max_run=max(runLen),mean_run=mean(runLen))

#check data for a particular bird
bird.names<-c("CI03H02_S84_L001")
bird.check<-lat.df%>%filter(name_in_vcf%in%bird.names)%>%
    group_by(name_in_vcf,ts,tagLifespan,tagDeployStart,recvDeployID,recvDeployLat,
             recvDeployLon,recvDeployName,detectyear,detectdoy)%>%
    mutate(ts=substr(ts,1,10))%>%mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%
    summarise(nDetect=n(),max_run=max(runLen),mean_run=mean(runLen))

ggplot(data = world) +
  geom_sf(fill="gray98")+
  geom_point(data=bird.check,
             aes(x=recvDeployLon,y=recvDeployLat,colour=detectdoy,
                             shape=as.factor(detectyear)),size=3)+
  coord_sf(xlim=c(-150,-40),ylim=c(10,60),expand=FALSE)


options(dplyr.print_min = Inf) 

#adjust so that start day is Aug 24 (day first bird tagged) not Jan 1
lat.df<-lat.df%>%
  mutate(tag.days.yr=ifelse(tagDays>365,tagDays-365,tagDays),
         yr=case_when(tagDays>365&tagDays<365*2~"year2",tagDays<366~"year1",tagDays>365*2~"year3"))%>%
  mutate(doy.aug24=case_when(ts>paste(detectyear,"08-24",sep="-")~difftime(ts,paste(detectyear,"08-24",sep="-"),units="days"),
                             ts<=paste(detectyear,"08-24",sep="-")~difftime(ts,paste((detectyear-1),"08-24",sep="-"),units="days")))%>%
  mutate(doy.aug24=round(as.numeric(gsub("days ","",doy.aug24))))%>%
  mutate(tag.days.yr=ifelse(tagDays>(365*2),tagDays-365*2,tag.days.yr))


print_maps=F
map_location="../../figures/filtered_birds_2023/check_231019/"

#print out a map of every included detection for each bird to manually check
if(print_maps==T){
  for(bird1 in unique(lat.df$name_in_vcf)){
    bird.lat<-lat.df%>%filter(name_in_vcf==bird1)%>%
      select(name_in_vcf,ts,tagDeployStart,recvDeployID,recvDeployLat,recvDeployLon,
             recvDeployName,detectyear,tagDays,detectdoy,tag.days.yr,doy.aug24,yr)%>%
      mutate(ts=substr(ts,1,10))%>%mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%
      distinct%>%arrange(tagDays)%>%
      filter(!is.na(tagDays))
    
    bird.title<-paste(bird1,"| towers:",length(unique(bird.lat$recvDeployName)),
                      "| days survived:",max(bird.lat$tagDays))
    vir.start<-min(bird.lat$doy.aug24)/365
    vir.end<-max(bird.lat$doy.aug24)/365
    
    png(paste(map_location,bird1,".png",sep=""),width=20,height=12,units="cm",res=300)
    print(ggplot(data = world) +
            geom_sf(fill="gray98")+
            geom_path(data=bird.lat,aes(x=recvDeployLon,y=recvDeployLat),size=0.8)+
            scale_colour_viridis(begin=vir.start,end=vir.end,
                                 name="day of year\nday one = Aug 24")+
            scale_shape_manual(values=unique(as.numeric(substr(bird.lat$yr,5,5))-1),
                               name="year")+
            geom_point(data=bird.lat,aes(x=recvDeployLon,y=recvDeployLat,
                                         colour=tag.days.yr,shape=yr),size=3,stroke=2)+
            theme_minimal()+ggtitle(bird.title)+
            coord_sf(xlim=c(-150,-40),ylim=c(2,60),expand=FALSE))
    
    dev.off() }}




#good example bird
bird.lat<-lat.df%>%filter(name_in_vcf=="CI04H02_S90_L001")%>%
  select(name_in_vcf,ts,tagDeployStart,recvDeployID,recvDeployLat,recvDeployLon,
         recvDeployName,detectyear,tagDays,detectdoy,tag.days.yr,doy.aug24)%>%
  mutate(ts=substr(ts,1,10))%>%mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%
  distinct%>%arrange(tagDays)%>%
  filter(!is.na(tagDays))

#png(paste("../figures/example_bird.png",sep=""),width=15,height=8,units="cm",res=300)
ggplot(data = world) +
  geom_sf(fill="gray95")+
  geom_path(data=bird.lat,aes(x=recvDeployLon,y=recvDeployLat),size=0.8)+
  scale_colour_viridis(name="days since tagging",option="C")+
  scale_shape_manual(values=unique(as.numeric(substr(bird.lat$yr,5,5))-1),name="year")+
  geom_point(data=bird.lat,aes(x=recvDeployLon,y=recvDeployLat,
                               colour=tagDays),size=3,stroke=2)+
  theme_minimal()+xlab("Longitude")+ylab("Latitude")+
  coord_sf(xlim=c(-150,-40),ylim=c(2,60),expand=FALSE)
dev.off()

#good example bird
birds.fall<-lat.df%>%filter(name_in_vcf%in%c("DH29H10","DI02H07","DI04H04"))%>%
  select(name_in_vcf,ts,tagDeployStart,recvDeployID,recvDeployLat,recvDeployLon,
         recvDeployName,detectyear,tagDays,detectdoy,tag.days.yr,doy.aug24)%>%
  mutate(ts=substr(ts,1,10))%>%mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%
  distinct%>%arrange(tagDays)%>%
  filter(!is.na(tagDays)&tagDays<180)

#png("../../figures/Fig1/Fig1B.png",width=18,height=14,units="cm",res=300)
ggplot() +
  geom_sf(data = world,fill="grey20",colour=NA)+
  geom_path(data=birds.fall,aes(x=recvDeployLon,y=recvDeployLat,colour=name_in_vcf),size=1.2)+
  scale_colour_viridis(discrete=T,begin=0.4,end=0.95)+
  #scale_colour_manual(values=c("lightsteelblue4","thistle4","honeydew4"))+
  geom_point(data=birds.fall,aes(x=recvDeployLon,y=recvDeployLat,
                               colour=name_in_vcf,shape=name_in_vcf),size=3,stroke=2)+
  scale_shape_manual(values=c(19,15,17))+
  xlab("Longitude")+ylab("Latitude")+
  theme(legend.position="none")+ylim(0,66)+xlim(-162,-55)
dev.off()

#png(paste("../figures/transect.png",sep=""),width=8,height=8,units="cm",res=300)
ggplot(data = world) +
  geom_sf(fill="gray95")+
  geom_point(data=lat.df%>%filter(recvProjID==unique(tagProjID))%>%
              select(recvDeployLat,recvDeployLon,recvDeployID)%>%distinct(),
            aes(x=recvDeployLon,y=recvDeployLat),size=3,shape=0,stroke=2,colour="darkslategray")+
  theme_minimal()+xlab("Longitude")+ylab("Latitude")+
  scale_x_continuous(breaks=seq(-126,-118,by=4))+
  scale_y_continuous(breaks=seq(48,52,by=2))+
  coord_sf(xlim=c(-126,-118),ylim=c(48,52),expand=FALSE)


#make legend for colour scale that starts on Aug 24
legend.df<-data.frame(plotDays=seq(as.Date("2020-08-24"),as.Date("2021-08-23"),by="days"),
           day_of_year=1:365)
legend.df$dayOfYear<-as.numeric(format(legend.df$plotDays, "%j"))

ggplot(legend.df,aes(y=dayOfYear,x=1,colour=day_of_year,shape=year))+
  geom_point(shape=15,size=10)+scale_colour_viridis()+
  scale_y_continuous(labels = function(x) format(as.Date(as.character(x), "%j"), "%d-%b"),
                     breaks=round(seq(min(legend.df$dayOfYear), max(legend.df$dayOfYear), by = 10),1))+
  theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+xlab("")+ylab("")

#plot location and time of all detections
lat.coord<-lat.df%>%
  select(name_in_vcf,ts,tagDeployStart,recvDeployID,recvDeployLat,recvDeployLon,
         recvDeployName,recvProjID,recvProjName,detectyear,detectdoy,tagDays,doy.aug24)%>%
  mutate(ts=substr(ts,1,10))%>%mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%distinct()
#write.csv(lat.coord,"../summarized_detections.csv",row.names=F)


ggplot(data = world) +
  geom_sf(fill="gray98")+
  geom_point(data=lat.coord,aes(x=recvDeployLon,y=recvDeployLat,colour=doy.aug24),shape=19,size=3)+
  theme_minimal()+scale_colour_viridis(name="day of year\nday one = Aug 24")+
  coord_sf(xlim=c(-150,-40),ylim=c(0,70),expand=FALSE)+
  xlab("Longitude")+ylab("Latitude")

###################################################################
#now sort out adult survival
#based on recapture of tagged birds
###################################################################

adults_recap<-rbind(adult_recap_20102022%>%
                      mutate(retrieved_archival=case_when(retrieved_archival=="y"~1,
                                                          retrieved_archival=="n"~0))%>%
                      rename(band=band_.)%>%
                      select(name_in_vcf,retrieved_archival,release_year,age_release,band),
                    meta_thrush%>%
                      #extract adults tagged in 2022
                      filter(age_release!="HY"&tag_type=="archival"&release_year==2022)%>%
                      mutate(retrieved_archival=0)%>%
                      select(name_in_vcf,retrieved_archival,release_year,age_release,band))

#fix typos then record birds retrieved in 2023 as recaptured
adults_recap<-adults_recap%>%
  filter(name_in_vcf!="AF14H02")%>% #remove duplicate individual (other ID DF27H08)
  mutate(band=if_else(name_in_vcf=="DF27H08","298108794",band))%>%
  mutate(retrieved_archival=if_else(band%in%adult_recap_2023$banding_.,1,retrieved_archival))%>%
  select(name_in_vcf,retrieved_archival) #only keep the two necessary variables

#add adult data to survival spreadsheet
HI_thrush<-left_join(HI_thrush,adults_recap)%>%
  filter(name_in_vcf!="AF14H02")


#write.csv(HI_thrush,"../thrush_survival_240106.csv")

