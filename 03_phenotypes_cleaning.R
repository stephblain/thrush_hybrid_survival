############################################
#Clean up raw phenotypes data
#
#Input: morphological, migratory, survival data to 2021 from Hannah,
##combined banding sheets, HIest estimates
#
#Output 1: combined phenotypes dataframe for radio tagged juvenile birds in HIest dataset
#If modifying to output adults as well, note that 2022 tagged adult migratory phenotypes not included here
#Output 2: metadata for adult and juvenile birds
############################################


#load packages
pkgs<-c("tidyverse","viridis","ggpubr")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

setwd("C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/raw/")

morph_all2021<-read.csv("phenotypes_morpho.csv") #morphological phenotypes 2010-2021
morph_geo2022<-read.csv("banding_2022_geos.csv") #combined banding sheet for 2022 geos
morph_rad2022<-read.csv("banding_2022_motus.csv") #combined banding sheet for 2022 radios

bear_rad2021<-read.csv("fall_bearing_radio_raw.csv") #raw fall bearings 2019-2021
time_rad2021<-read.csv("radio_timing_raw.csv") #raw fall timing 2019-2021
migrate_rad2022<-read.csv("all_new_phenotypes_fall_2022.csv") #fall timing and bearings 2022

meta_all2021<-read.csv("survival_data_motus_by_latitude_updated.csv") #metadata 2010-2021
refs_all2021<-read.csv("order_stitch.csv") #references and name_in_vcf
sexChrom_all2022<-read.table("sexChrom_depths.txt",header=T) #read depths of sex chroms
sexChrom_unclear<-read.csv("sex_unclear_kd.csv",row.names="X") #sexes checked in lab

meta_geo2022<-morph_geo2022 #set aside for metadata
meta_rad2022<-morph_rad2022 #set aside for metadata

HI_thrush<-read.csv("hiest.stitch.ldr01.hwe.s1.fst94.noZ.adj98.csv") #S&H estimates

############################################
#morphological phenotypes
############################################

morph_all2021<-morph_all2021%>%rename_with(tolower)%>%
  rename(band=band_.,kipps=ps_kipps)%>%
  filter(name_in_vcf%in%HI_thrush$ind)%>%
  select(-order,-sex_binary)%>%filter(reference!="")

morph_geo2022<-morph_geo2022%>%rename_with(tolower)%>%
  rename(reference=reference..,band=band..,release_site=site,
         age_release=age,
         kipps=p...s,geo=geo..,nest=nest..)%>%
  mutate(name_in_vcf=reference)%>%
  filter(name_in_vcf%in%HI_thrush$ind)%>%
  select(colnames(morph_all2021))

morph_rad2022<-morph_rad2022%>%rename_with(tolower)%>%
  rename(reference=reference..,band=band..,release_site=site,
         age_release=age,
         kipps=p...s,geo=geo..,nest=nest..)%>%
  mutate(name_in_vcf=reference)%>%
  filter(name_in_vcf%in%HI_thrush$ind)%>%
  select(colnames(morph_all2021))

morph_all<-rbind(morph_all2021,morph_geo2022,morph_rad2022)

if(length(unique(morph_all$reference))==nrow(morph_all)){
 rm(morph_geo2022,morph_rad2022,morph_all2021)}else{print("check row binding")}


#correct errors and remove biologically improbable outliers
#errors identified by comparison to banding sheets
morph_all$p9[morph_all$name_in_vcf=="AF27H04"]<-74 #change from 7.4
morph_all$tail.length[morph_all$name_in_vcf=="AF20H01"]<-69
morph_all$tarsus.length[morph_all$name_in_vcf%in%c("AH26H03_S158_L001","BI06H03")]<-NA
morph_all$distal[morph_all$name_in_vcf%in%c("CI06H04_S107_L001","CG01H02_S31_L001")]<-NA
morph_all$kipps[morph_all$name_in_vcf=="BH31H12"]<-28.7
morph_all$tail.length[morph_all$name_in_vcf=="DI04H03"]<-72 #two values switched in entering data
morph_all$tarsus.length[morph_all$name_in_vcf=="DI04H03"]<-29
morph_all$p9[morph_all$name_in_vcf=="DI01H15"]<-NA #value entered 84 - wrong based on p7, p8, weird banding entry



############################################
#migratory phenotypes
############################################

bear_rad2021<-bear_rad2021%>%rename_with(tolower)%>%
  rename(motustagid=id)
time_rad2021<-time_rad2021%>%rename_with(tolower)
migrate_rad2021<-left_join(bear_rad2021,time_rad2021)
if(nrow(migrate_rad2021)==nrow(bear_rad2021)){
  rm(time_rad2021,bear_rad2021)}else{print("check left join")}
migrate_rad2022<-migrate_rad2022%>%
  rename(motustagid=id)%>%
  mutate(reference=substr(name_in_vcf,1,7))%>%
  select(-grep("qnorm",colnames(migrate_rad2022)))
migrate_rad<-rbind(migrate_rad2021%>%select(colnames(migrate_rad2022)),
                   migrate_rad2022)


############################################
#collect metadata
############################################



##no name in vcf column in morph all
#use banding sheets

meta_all2021<-left_join(meta_all2021,refs_all2021%>%select(name_in_vcf,reference))
meta_all2021<-meta_all2021%>%as_tibble()%>%
  rename_with(tolower)%>%
  rename(band=band_.)%>%
  filter(name_in_vcf%in%HI_thrush$ind)%>%
  select("reference","name_in_vcf","release_site","release_year","tag_type",
         "sex_binary","age_release","release_gps.n",
         "release_gps.w","motustagid","band")

meta_geo2022<-meta_geo2022%>%as_tibble()%>%
  rename_with(tolower)%>%
  rename(reference=reference..,release_site=site,sex_binary=sex,age_release=age,
         release_gps.n=gps.n,release_gps.w=gps.w,band=band..)%>%
  mutate(sex_binary=case_when(sex_binary=="M"~"1",sex_binary=="F"~"0"))%>%
  mutate(release_year=2022,tag_type="archival",motustagid=NA,name_in_vcf=reference)%>%
  filter(name_in_vcf%in%HI_thrush$ind)%>%
  mutate(age_release=gsub("check pictures","AHY",age_release))%>%
  mutate(age_release=gsub("\\?","",age_release))%>%
  select(colnames(meta_all2021))

meta_rad2022<-meta_rad2022%>%as_tibble()%>%
  rename_with(tolower)%>%
  rename(reference=reference..,release_site=site,sex_binary=sex,age_release=age,
         release_gps.n=gps.n,release_gps.w=gps.w,band=band..)%>%
  mutate(sex_binary=case_when(sex_binary=="M"~"1",sex_binary=="F"~"0"))%>%
  mutate(release_year=2022,tag_type="radio",name_in_vcf=reference)%>%
  filter(name_in_vcf%in%HI_thrush$ind)%>%
  select(colnames(meta_all2021))

#get sex estimates for 2022 birds

sexChrom_all2022<-sexChrom_all2022%>%
  left_join(rbind(meta_rad2022%>%select(name_in_vcf,tag_type),
                  meta_geo2022%>%select(name_in_vcf,tag_type)))%>%
  drop_na()

sexChrom_all2022.l<-sexChrom_all2022%>%
  mutate(scaf7=scaf7/scaf3,scaf25=scaf25/scaf3)%>%
  select(-scaf3)%>%
  #mutate(sex_maybe=case_when(scaf7<0.8~"female",scaf7>0.8&scaf7<1.0~"unclear",
  #                           scaf7>1.0~"male"))%>%
  mutate(sex_maybe=case_when(scaf7<1~"female",scaf7>1.0~"male"))%>%
  mutate(sex_maybe=if_else(tag_type=="archival","male",sex_maybe))%>%
  pivot_longer(cols=c(scaf7,scaf25),names_to="scaffold",
               values_to="fold_diff")

#replace unknown sexes with pcr-inferred sex
sexChrom_all2022.l<-sexChrom_all2022.l%>%
  left_join(sexChrom_unclear%>%select(name_in_vcf,sex)%>%rename(sex_PCR=sex))%>%
  mutate(sex_PCR=case_when(sex_PCR=="F"~"female",sex_PCR=="M"~"male",TRUE~NA))%>%
  mutate(sex_maybe=if_else(!is.na(sex_PCR),sex_PCR,sex_maybe))%>%
  select(-sex_PCR)%>%rename(sex=sex_maybe)


ggplot(sexChrom_all2022.l,aes(x=fct_reorder(name_in_vcf,fold_diff),
                            y=fold_diff,colour=sex,group=name_in_vcf,
                            shape=scaffold))+
  geom_hline(yintercept=c(0.8,1))+
  geom_line(colour="grey30")+geom_point()+
  scale_colour_manual(values=c("coral3","cyan4","grey30"))+
  ylab("fold diff read depth\nsex chrom / autosome")+xlab("reference")+
  theme(axis.text.x=element_text(angle = 90,vjust = 1,hjust=1))



ggplot(sexChrom_all2022.l,aes(x=fold_diff,fill=sex))+
  geom_histogram()+
  scale_fill_manual(values=c("coral3","cyan4","grey30"))+
  xlab("fold diff read depth\nsex chrom / autosome")+
  facet_grid(rows=vars(scaffold))+
  theme_minimal()


sexChrom_all2022.l<-sexChrom_all2022.l%>%
  filter(scaffold=="scaf7"&name_in_vcf%in%meta_rad2022$name_in_vcf)%>%
  mutate(sex_binary=case_when(sex=="male"~"1",sex=="female"~"0"))

if(all(sexChrom_all2022.l$name_in_vcf==meta_rad2022$name_in_vcf)){
  meta_rad2022$sex_binary<-sexChrom_all2022.l$sex_binary}else{
    print("check name_in_vcf columns")}


meta_all<-rbind(meta_all2021,meta_geo2022,meta_rad2022)


#individuals with no meta data are the montana birds - no pheno data for them anyway
#HI_thrush[!HI_thrush$ind%in%meta_all$name_in_vcf,]

if(length(unique(meta_all$reference))==nrow(meta_all)){
  rm(meta_geo2022,meta_rad2022,meta_all2021,refs_all2021)}else{
    print("check row binding")}


############################################
#combine all
############################################

HI_thrush<-HI_thrush%>%
  rename(ancestry=S,heterozygosity=H,name_in_vcf=ind)%>%
  select(ancestry,heterozygosity,name_in_vcf)

#missing phenotypes are for captive birds only
meta_all%>%
  filter(name_in_vcf%in%HI_thrush[!HI_thrush$name_in_vcf%in%morph_all$name_in_vcf,"name_in_vcf"])%>%
  pull(tag_type)%>%unique()

#only 1 radio tagged individual that is after second year
meta_all%>%
  group_by(tag_type,age_release)%>%
  summarize(count=n())



#keep only radio-tagged juvies from Pemberton
phenotypes<-HI_thrush%>%
  left_join(meta_all)%>%left_join(morph_all)%>%left_join(migrate_rad)%>%
  filter(tag_type=="radio"&age_release=="HY"&release_site=="Pemberton")%>%
  mutate(fat.score=as.numeric(fat.score),p9=as.numeric(p9))%>%
  #remove birds translocated from whistler
  filter(!name_in_vcf%in%c("BH30H01","BH29H01","BH29H03","BH29H04","BH29H05","BH29H07"))

#write to file
# write.csv(phenotypes,row.names=F,
#           "C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/thrush_phenotypes_HY_Pemberton.240104.csv")


############################################
#organize and print metaData
############################################

meta_all<-left_join(HI_thrush,meta_all)

#add info for Montana birds
meta_all<-rbind(meta_all%>%filter(!startsWith(name_in_vcf,"2")),
                meta_all%>%filter(startsWith(name_in_vcf,"2"))%>%
                  mutate(reference=name_in_vcf,release_site="Montana",tag_type="no_tag"))
all(!is.na(meta_all$tag_type))

#write to file
# write.csv(meta_all,row.names=F,
#           "C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/thrush_meta_240104.csv")
