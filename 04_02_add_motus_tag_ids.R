##################################################################################################
# script 2
##################################################################################################

# 1) Load libraries

library(dplyr)

##################################################################################################

# 2) Load meta data

# get information of motus tag deployment id
motus_id<-read.csv("C:/Users/hcjusten/Dropbox/PhD/Field/2022/MotusTagsTemplate-280.csv")
motus_id<-data.frame(motus_id$tagID,motus_id$tagProjectID,motus_id$mfgID,motus_id$tagDeployID,motus_id$utcYearStart)
names(motus_id)=c("motusTagID","tagProjectID","mfgID","tagDeployID","utcYearStart")

# get info of the individuals included in the survival dataset (connection between reference and name in vcf)
meta_info<-read.csv("C:/Users/hcjusten/Dropbox/PhD/Thesis/Chapter_2/Analysis/GWAS/meta_data/meta_stitch.csv")
meta_info<-data.frame(meta_info$order,meta_info$reference,meta_info$name_in_vcf)
names(meta_info)=c("order","reference","name_in_vcf")

# get info from banding sheets (conncetion between reference and which tag was attached to which bird)

# data from 2019
band19<-read.csv("C:/Users/hcjusten/Documents/delmore_lab/SWTH/banding_sheets/banding_2019_motus.csv")
band19<-data.frame(band19$Reference..,band19$Site,band19$mfgID,band19$motusTagID,band19$tagDeployID)
names(band19)=c("reference","release_site","mfgID","motusTagID","tagDeployID")

# data from 2020
band20<-read.csv("C:/Users/hcjusten/Documents/delmore_lab/SWTH/banding_sheets/banding_2020_motus.csv")
band20<-data.frame(band20$Reference..,band20$Site,band20$mfgID,band20$motusTagID,band20$tagDeployID)
names(band20)=c("reference","release_site","mfgID","motusTagID","tagDeployID")

# data from 2021
band21<-read.csv("C:/Users/hcjusten/Documents/delmore_lab/SWTH/banding_sheets/banding_2021_motus.csv")
band21<-data.frame(band21$Reference..,band21$Site,band21$mfgID,band21$motusTagID,band21$tagDeployID)
names(band21)=c("reference","release_site","mfgID","motusTagID","tagDeployID")

##################################################################################################

# 3) Merge motus data with meta data

motus_id19<-subset(motus_id,motus_id$utcYearStart==2019)
motus_id_sub19<-join(motus_id19,band19,by="motusTagID",match="first")

motus_id20<-subset(motus_id,motus_id$utcYearStart==2020)
motus_id_sub20<-join(motus_id20,band20,by="motusTagID",match="first")

motus_id21<-subset(motus_id,motus_id$utcYearStart==2021)
motus_id_sub21<-join(motus_id21,band21,by="motusTagID",match="first")

motus_id_final<-rbind(motus_id_sub19,motus_id_sub20,motus_id_sub21)

# add vcf reference to data set

meta_data<-join(motus_id_final,meta_info,by="reference",match="first")

##################################################################################################

df.alltags_sub_merged<-join(df.alltags.sub1,meta_data,by="motusTagID",match="first")

write.csv(df.alltags_sub_merged,"filtered_data_Nov_2022.csv",row.names=F)

##################################################################################################
# end of script
##################################################################################################