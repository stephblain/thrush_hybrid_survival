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

map_base<-ggplot() +
  geom_sf(data=world,fill="grey10",colour=NA)+
  geom_sf(data = coastal[1:2,1],fill=sp.cols[1],colour=sp.cols[1])+
  geom_sf(data = inland[1:2,1],fill=sp.cols[2],colour=sp.cols[2])


receivers<-receivers%>%
  mutate(dtEnd=as.Date(substr(dtEnd,1,10)))%>%
  filter(dtEnd>"2019-08-01"|deploymentStatus=="active")


#png(filename="../figures/Fig1.png",units="px")
A.main<-map_base+
  geom_hline(yintercept=c(25,40),colour="grey50",linewidth=1)+
  geom_point(data=receivers,aes(y=latitude,x=longitude),colour="maroon",size=0.7)+
  geom_polygon(aes(x=c(-125,-125,-119,-119),y=c(48.5,51.5,51.5,48.5)),
               fill=NA,colour="black",size=0.8)+
  xlab("Longitude")+ylab("Latitude")+
  ylim(0,66)+xlim(-162,-55)+
  theme(panel.background=element_rect(colour="black",size=1))

A.inset<-map_base+
  geom_point(data=lat.coord%>%filter(recvProjName=="thrushes")%>%
               select(recvDeployLat,recvDeployLon,recvDeployID)%>%distinct(),
             aes(x=recvDeployLon,y=recvDeployLat),size=2,shape=19,
             colour="maroon")+
  coord_sf(xlim=c(-125,-119),ylim=c(48.5,51.5),expand=FALSE)+
  theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
        panel.border=element_rect(colour="black",fill=NA))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))


B.main<-map_base+
  geom_line(data=lat.coord%>%filter(doy.aug24<183)%>%
              mutate(groupID=paste(name_in_vcf,detectyear)),
            aes(y=recvDeployLat,x=recvDeployLon,group=groupID),colour="maroon",
            size=0.5)+
  geom_polygon(aes(x=c(-125,-125,-119,-119),y=c(48.5,51.5,51.5,48.5)),
               fill=NA,colour="black",size=0.8)+
  xlab("Longitude")+ylab("Latitude")+
  ylim(0,66)+xlim(-162,-55)+
  theme(panel.background=element_rect(colour="black",size=1))

B.inset<-map_base+
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



#make tri plots too
#read in file with HIest estimates, individual names, and survival data
HI_thrush<-read.csv("C:/Users/Steph/GitHub_data/hiest/thrush_survival_231130.csv"); HI_thrush$X<-NULL

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

HI_adult<-HI_thrush%>%filter(age_release%in%c("AHY","SY","ASY"))%>%
  mutate(retrieved_archival=as.numeric(retrieved_archival))%>%
  mutate(age_class="adult",survival=retrieved_archival)%>%
  filter(!is.na(retrieved_archival))
HI_juvie<-HI_thrush%>%filter(age_release=="HY"&pop_type=="hybrid"&release_site=="Pemberton")%>%
  mutate(age_class="juvenile",survival=t5_spring40)
HI_both<-rbind(HI_adult,HI_juvie)%>%
  mutate(pop_type=if_else(pop_type=="hybrid","1hybrid",pop_type))%>%
  arrange(pop_type)


tri<-ggplot(data=HI_both,aes(x=ancestry,y=heterozygosity))+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="grey30",size=0.8)+
  #geom_line(aes(x =ancestry, y = 0,colour=ancestry),size=3)+
  geom_point(aes(colour=as.factor(pop_type)),size=0.7,stroke=1,shape=19)+
  scale_colour_manual(values=c("grey10",sp.cols[1],sp.cols[2]))+
  #scale_colour_viridis(begin=0.4,end=0.95,guide="none")+
  facet_grid(cols=vars(age_class))+
  theme(strip.background = element_rect(colour="grey10", fill="white", 
                                        size=1, linetype="solid"),
        strip.text.x = element_text(size=11))+
  xlab("Ancestry")+ylab("Heterozygosity")+coord_equal()+
  theme(panel.background=element_rect(colour="black",size=1),legend.position="none")

tri

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