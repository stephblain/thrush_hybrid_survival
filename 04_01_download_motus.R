##################################################################################################
#
# script to prepare survival data
#
##################################################################################################

# 1) Load libraries

library(remotes)
library(motus)
library(motusData)
library(maps)
library(tidyverse)
library(rworldmap)
library(ggmap)
library(lubridate)
library(curl)
library(DBI)
library(RSQLite)
library(dbplyr)

##################################################################################################
##################################################################################################

# 2) Load data

# set working directory, data will be stored here

setwd("C:/Users/hcjusten/Dropbox/PhD/Thesis/Chapter_2/Analysis/Motus/data/july2_2021")


Sys.setenv(TZ = "UTC")

# working with sample data project

proj.num <- 280

# when downloading data for the first time: you will be asked for your motus name and password

sql.motus <- tagme(projRecv = proj.num, new = F, update = TRUE)

##################################################################################################
##################################################################################################
# 3) filter data

#motus data filtered:
tbl_motusFilter <- filterByActivity(sql.motus, return = "all")

dat<-subset(tbl_motusFilter, tbl_motusFilter$probability==1)
#dat<-data.frame(subset(tbl_motusFilter,tbl_motusFilter$runLen>3))
 write.csv(dat, file="dat.csv",row.names=F)

df.alltags.sub <- tbl_motusFilter %>% 
  filter(probability == 1) %>%
  collect() %>%
  as.data.frame() %>%
  mutate(ts = as_datetime(ts),  # work with dates AFTER transforming to flat file
         tagDeployStart = as_datetime(tagDeployStart),
         tagDeployEnd = as_datetime(tagDeployEnd))

df.alltags.sub_year <- df.alltags.sub %>%
  mutate(ts = as_datetime(ts, tz = "UTC"), # convert ts to POSIXct format
         year = year(ts), # extract year from ts
         doy = yday(ts)) %>% # extract numeric day of year from ts
  filter(!is.na(recvDeployLat))

##################################################################################################
##################################################################################################

# additional filtering:

# subset to tags registered to our motus project
df.alltags.sub.merged1<-subset(df.alltags.sub.merged,df.alltags.sub.merged$tagProjID==280)

# remove detections from towers taht frequenly show wrong detections and are out of the range where we expect the birds
df.alltags.sub.merged1 <- subset(df.alltags.sub.merged1, df.alltags.sub.merged1$recvDeployName != "Tadoussac - Andr?")
df.alltags.sub.merged1 <- subset(df.alltags.sub.merged1, df.alltags.sub.merged1$recvDeployName != "Lockoff")
df.alltags.sub.merged1 <- subset(df.alltags.sub.merged1, df.alltags.sub.merged1$recvDeployName != "Triton2")

# remove detections that are more east than 70°W
df.alltags.sub1<- subset(df.alltags.sub.merged1, df.alltags.sub.merged1$recvDeployLon < -70)

##################################################################################################
##################################################################################################

write.csv(df.alltags.sub1,"filtered_motus_data.csv",row.names=F)

##################################################################################################
# end of script
##################################################################################################
