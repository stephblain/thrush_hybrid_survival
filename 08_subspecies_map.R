#load packages
pkgs<-c("tidyverse","viridis","ggpubr","sf","rnaturalearth","ggspatial","cowplot")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

setwd("C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/")

lat.coord<-read.csv("summarized_detections.csv")
coastal<-read_sf(dsn="C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/figures/Fig1/Catharus_shpFiles",
                   layer="Catharus_ustulatus")
inland<-read_sf(dsn="C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/figures/Fig1/Catharus_shpFiles",
                 layer="Catharus_swainsoni")
receivers<-read.csv("receiver-deployments.csv")

#load map data
world <- ne_countries(scale = "medium", returnclass = "sf")

sp.cols<-viridis(2,begin=0.4,end=0.95)

# map_main1<-ggplot() +
#   geom_sf(data=world,fill="gray95")+
#   geom_sf(data = coastal[1:2,1],fill="lightsteelblue3")+
#   geom_sf(data = inland[1:2,1],fill="honeydew3")+
#   ylim(0,66)+xlim(-162,-55)+
#   geom_polygon(aes(x=c(-125,-125,-119,-119),y=c(48.5,51.5,51.5,48.5)),
#                fill=NA,colour="black")+
#   xlab("Longitude")+ylab("Latitude")
# 
# map_main<-ggplot() +
#   geom_sf(data=world,fill="grey20",colour=NA)+
#   geom_sf(data = coastal[1:2,1],fill=sp.cols[1],colour=sp.cols[1])+
#   geom_sf(data = inland[1:2,1],fill=sp.cols[2],colour=sp.cols[2])+
#   ylim(0,66)+xlim(-162,-55)+
#   geom_polygon(aes(x=c(-125,-125,-119,-119),y=c(48.5,51.5,51.5,48.5)),
#                fill=NA,colour="black",size=0.8)+
#   xlab("Longitude")+ylab("Latitude")
# 
# 
# map_inset1<-ggplot() +
#   geom_sf(data=world,fill="gray95")+
#   geom_sf(data = coastal[1:2,1],fill="lightsteelblue3")+
#   geom_sf(data = inland[1:2,1],fill="honeydew3")+
#   geom_point(data=lat.coord%>%filter(recvProjName=="thrushes")%>%
#                select(recvDeployLat,recvDeployLon,recvDeployID)%>%distinct(),
#              aes(x=recvDeployLon,y=recvDeployLat),size=3,shape=19,
#              colour="tomato3")+
#   coord_sf(xlim=c(-125,-119),ylim=c(48.5,51.5),expand=FALSE)+
#   theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
#         panel.border=element_rect(colour="black",fill=NA))+
#   theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))



map_base<-ggplot() +
  geom_sf(data=world,fill="grey10",colour=NA)+
  geom_sf(data = coastal[1:2,1],fill=sp.cols[1],colour=sp.cols[1])+
  geom_sf(data = inland[1:2,1],fill=sp.cols[2],colour=sp.cols[2])

# map_main<-map_base+xlab("Longitude")+ylab("Latitude")+
#   ylim(0,66)+xlim(-162,-55)+
#   geom_polygon(aes(x=c(-125,-125,-119,-119),y=c(48.5,51.5,51.5,48.5)),
#                fill=NA,colour="black",size=0.8)
# 
# 
# map_inset<-map_base+
#   geom_point(data=lat.coord%>%filter(recvProjName=="thrushes")%>%
#                select(recvDeployLat,recvDeployLon,recvDeployID)%>%distinct(),
#              aes(x=recvDeployLon,y=recvDeployLat),size=3,shape=19,
#              colour="maroon4")+
#   coord_sf(xlim=c(-125,-119),ylim=c(48.5,51.5),expand=FALSE)+
#   theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
#         panel.border=element_rect(colour="black",fill=NA))+
#   theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))


# #png("../figures/Fig1/Fig1A.png",width=18,height=14,units="cm",res=300)
# ggdraw() +
#   draw_plot(map_main) +
#   draw_plot(map_inset, x = 0.115, y = .15, width = .3, height = .3)
# dev.off()

receivers<-receivers%>%
  mutate(dtEnd=as.Date(substr(dtEnd,1,10)))%>%
  filter(dtEnd>"2019-08-01"|deploymentStatus=="active")


#png(filename="../figures/Fig1.png",units="px")
A.main<-map_base+
  geom_hline(yintercept=c(25,40),colour="grey50",size=1)+
  geom_point(data=receivers,aes(y=latitude,x=longitude),colour="maroon",size=0.7)+
  geom_polygon(aes(x=c(-125,-125,-119,-119),y=c(48.5,51.5,51.5,48.5)),
               fill=NA,colour="black",size=0.8)+
  xlab("Longitude")+ylab("Latitude")+
  ylim(0,66)+xlim(-162,-55)

A.inset<-map_base+
  geom_point(data=lat.coord%>%filter(recvProjName=="thrushes")%>%
               select(recvDeployLat,recvDeployLon,recvDeployID)%>%distinct(),
             aes(x=recvDeployLon,y=recvDeployLat),size=2,shape=19,
             colour="maroon")+
  coord_sf(xlim=c(-125,-119),ylim=c(48.5,51.5),expand=FALSE)+
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        panel.border=element_rect(colour="black",fill=NA))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))


dev.off()
B.main<-map_base+
  geom_line(data=lat.coord%>%filter(doy.aug24<183)%>%
              mutate(groupID=paste(name_in_vcf,detectyear)),
            aes(y=recvDeployLat,x=recvDeployLon,group=groupID),colour="maroon",
            size=0.5)+
  geom_polygon(aes(x=c(-125,-125,-119,-119),y=c(48.5,51.5,51.5,48.5)),
               fill=NA,colour="black",size=0.8)+
  xlab("Longitude")+ylab("Latitude")+
  ylim(0,66)+xlim(-162,-55)

B.inset<-map_base+
  # geom_point(data=lat.coord%>%filter(recvProjName=="thrushes")%>%
  #              select(recvDeployLat,recvDeployLon,recvDeployID)%>%distinct(),
  #            aes(x=recvDeployLon,y=recvDeployLat),size=3,shape=19,
  #            colour="maroon")+
  geom_line(data=lat.coord%>%filter(doy.aug24<183)%>%
              mutate(groupID=paste(name_in_vcf,detectyear)),
            aes(y=recvDeployLat,x=recvDeployLon,group=groupID),colour="maroon",
            size=0.5)+
  coord_sf(xlim=c(-125,-119),ylim=c(48.5,51.5),expand=FALSE)+
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        panel.border=element_rect(colour="black",fill=NA))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))



#png("../figures/Fig1/Fig1.png",width=18,height=7,units="cm",res=300)
ggarrange(ggdraw() +
  draw_plot(A.main) +
  draw_plot(A.inset, x = 0.18, y = .2, width = .3, height = .3),
  ggdraw() +
    draw_plot(B.main) +
    draw_plot(B.inset, x = 0.18, y = .2, width = .3, height = .3),
  labels=LETTERS)
#dev.off()

ggarrange(A,B,labels=LETTERS)

#make tri plots too
#read in file with HIest estimates, individual names, and survival data
HI_thrush<-read.csv("thrush_survival_230712.csv"); HI_thrush$X<-NULL

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

HI_adult<-HI_thrush%>%filter(age_release%in%c("AHY","SY","ASY")&pop_type=="hybrid")%>%
  mutate(retrieved_archival=as.numeric(retrieved_archival))%>%
  mutate(age_class="adult",survival=retrieved_archival)%>%
  filter(!is.na(retrieved_archival))
HI_juvie<-HI_thrush%>%filter(age_release=="HY"&pop_type=="hybrid"&release_site=="Pemberton")%>%
  mutate(age_class="juvenile",survival=t5_spring40)
HI_both<-rbind(HI_adult,HI_juvie)
HI_both$survival
tri<-ggplot(data=HI_both,aes(x=ancestry,y=heterozygosity))+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey50",size=2)+
  geom_point(aes(shape=as.factor(survival)),colour="darkslategray",size=2,stroke=1)+
  scale_shape_manual(values=c(1,19),guide="none")+
  facet_grid(cols=vars(age_class))+
  theme(strip.background = element_rect(colour="grey50", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size=11))+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()

viridis(begin=0.67,end=0.68,n=1)
(0.95-0.4)/2

tri<-ggplot(data=HI_both,aes(x=ancestry,y=heterozygosity))+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey10",size=2)+
  geom_point(aes(shape=as.factor(survival),colour=ancestry),size=2,stroke=1)+
  scale_shape_manual(values=c(1,19),guide="none")+
  scale_colour_viridis(begin=0.4,end=0.95,guide="none")+
  facet_grid(cols=vars(age_class))+
  theme(strip.background = element_rect(colour="grey10", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size=11))+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()


tri<-ggplot(data=HI_both,aes(x=ancestry,y=heterozygosity))+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey70",size=3)+
  #geom_line(aes(x =ancestry, y = 0,colour=ancestry),size=3)+
  geom_point(aes(shape=as.factor(survival)),colour="grey10",size=2,stroke=1)+
  scale_shape_manual(values=c(1,19),guide="none")+
  #scale_colour_viridis(begin=0.4,end=0.95,guide="none")+
  facet_grid(cols=vars(age_class))+
  theme(strip.background = element_rect(colour="grey10", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size=11))+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()
blank<-ggplot()+geom_blank()+
  theme(axis.line=element_line(color="white"))

#png("../figures/Fig1/Fig1.png",width=18,height=14,units="cm",res=300)
#pdf("../figures/Fig1/Fig1.pdf",width=7.09,height=5.51)
gg1<-ggarrange(ggdraw() +
            draw_plot(A.main) +
            draw_plot(A.inset, x = 0.18, y = .2, width = .3, height = .3),
          ggdraw() +
            draw_plot(B.main) +
            draw_plot(B.inset, x = 0.18, y = .2, width = .3, height = .3),
          labels=c("A","B"))
gg2<-ggarrange(tri,blank,widths=c(0.7,0.3),labels=c("C","D"))
ggarrange(gg1,gg2,nrow=2)

dev.off()

tracks<-map_base+
  geom_line(data=lat.coord%>%filter(doy.aug24<183)%>%
              mutate(groupID=paste(name_in_vcf,detectyear)),
            aes(y=recvDeployLat,x=recvDeployLon,group=groupID),colour="maroon",
            size=0.5)+
  xlab("Longitude")+ylab("Latitude")+
  ylim(0,66)+xlim(-162,-55)

#png("../figures/Fig1/presentation_maps.png",width=18,height=7,units="cm",res=300)
ggarrange(map_base+
            xlab("Longitude")+ylab("Latitude")+
            ylim(0,66)+xlim(-162,-55),
          tracks)
#dev.off()



p1.tri<-ggplot()+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey10",size=2)+
  geom_point(data=HI_adult,
             aes(x=ancestry,y=heterozygosity,
                 shape=as.factor(retrieved_archival)),colour="darkslategrey",size=2,stroke=1)+
  scale_shape_manual(values=c(1,19),guide="none")+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()


p2.tri<-ggplot()+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey10",size=2)+
  geom_point(data=HI_juvie,
             aes(x=ancestry,y=heterozygosity,
                 shape=as.factor(t5_spring40)),colour="darkslategrey",size=2,stroke=1)+
  scale_shape_manual(values=c(1,19),guide="none")+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()

ggarrange(p1.tri,p2.tri,labels=LETTERS)
