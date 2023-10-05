############################################
#Build a network of juvenile phenotypes
############################################

#load packages
pkgs<-c("tidyverse","viridis","ggpubr","qgraph")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

setwd("C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/")


phenotypes<-read.csv("thrush_phenotypes_HY_Pemberton.csv")

colSums(!is.na(phenotypes))

phenoInput<-phenotypes%>%
  select(doy_fall_r1,bearing_fall_1,kipps,
         tail.length,tarsus.length,wing.cord)%>%
  drop_na()

migrateInput<-phenotypes%>%
  select(doy_fall_r1,doy_fall_r2,doy_fall_r4,
         bearing_fall_1,bearing_fall_2,bearing_fall_4,
         distal,p9, #delet fat score because minimal variation
         p10,kipps,tail.length,tarsus.length,wing.cord)%>%
  drop_na()

hiInput<-phenotypes%>%
  select(distal,doy_fall_r1,bearing_fall_1,p9,
         p10,kipps,tail.length,tarsus.length,wing.cord,
         ancestry)%>%
  drop_na()

#----------------------------------
## Functions from Wilkins et al. 2015
#A. Function to get confidence intervals for bootstraps of correlation matrix

boot.cor<-function(DATA,shps,iterations)#DATA is correlation matrix, shps are modality designations, iterations is number of bootstraps
{
  t=iterations # number of iterations for bootstrapping; 
  pb <- txtProgressBar(min=1, max=t-1,style=3) #define progress bar
  cor.orig=cor(DATA,method="spearman",use="pairwise.complete.obs")
  diag(cor.orig)=0
  boot.cor=array(dim=c(nrow(cor.orig),ncol(cor.orig),t))
  cors=vector(length=t)
  #Loop to create t=iterations datasets composed of k trait observations from n individuals resampled with replacement from our original dataset
  for (i in 1:t)
  {
    setTxtProgressBar(pb,i)#update progress bar
    s=sample(1:nrow(DATA),nrow(DATA),replace=T) #sample rows with replacement
    new.dat=DATA[s,] #use the row numbers from above to construct bootstrapped sample
    m=cor(new.dat,method="spearman",use="pairwise.complete.obs") #new correlation matrix
    boot.cor[,,i]=m
    cors[i]=cor.test(m[which(upper.tri(m)==TRUE)],cor.orig[which(upper.tri(cor.orig)==TRUE)])$estimate
    #       assort.boot.orig[i]=assortment.discrete(abs(m),shps)$r #calculate assortativity. Use absolute values.
    #       assort.boot.filt[i]=assortment.discrete(abs(m.filt),shps)$r #calculate assortativity for filtered bootstrap data matrix
  }
  
  avg.cor=apply(boot.cor,c(1,2),mean,na.rm=T) #calculate average correlation matrix from 1000 resampling iterations
  #calculate bootstrap confidence intervals for each edge in the original dataset
  lower.cor=apply(boot.cor, c(1,2), quantile, probs=c(0.025))
  upper.cor=apply(boot.cor, c(1,2), quantile, probs=c(0.975))
  return(list(orig.cor=cor.orig,avg.cor=avg.cor,lower.cor=lower.cor,upper.cor=upper.cor))
}
######### End bootstrapping function


######
# B. Function to plot bootstrap correlation CIs vs empirical correlations, and extract robust edges

#Analyze output of boot.cor function, generating assorativity and redundancy metrics
#Also generates figure showing bootstrap correlations against empirical correlations

analyze.boot<-function(boot.cor.output)
{
  par(cex.lab=1.8,mar=c(5,5,3,3),cex.axis=1.4,font.lab=2,mgp=c(3.5,1,0))  
  D<-boot.cor.output
  emp.edges=D$orig.cor[which(upper.tri(D$orig.cor)==TRUE)]
  boot.avgs=D$avg.cor[which(upper.tri(D$avg.cor)==TRUE)]
  boot.lower=D$lower.cor[which(upper.tri(D$lower.cor)==TRUE)]
  boot.upper=D$upper.cor[which(upper.tri(D$upper.cor)==TRUE)]
  ci.0=c("black","red")[((boot.lower>0&boot.upper>0)|(boot.lower<0&boot.upper<0))+1]
  plot(emp.edges,boot.avgs, pch=19, ylim=c(-1,1), xlim=c(-1,1), xlab="Empirical Correlations", ylab="Bootstrap Correlations",col=ci.0,las=1)
  for (i in 1:length(emp.edges)){
    lines(c(emp.edges[i], emp.edges[i]), c(boot.lower[i],boot.upper[i]),col=ci.0[i])
  }
  min.threshold<-
    abline(v=0.3,lty=2)
  abline(v=-0.3,lty=2)
  abline(h=0, lty=1)
  
  # find threshold where bootstrap CIs don't overlap 0.
  abs(boot.lower)>0
  mini=maxi=0
  for(i in 1: length(emp.edges))
  {
    if(emp.edges[i]<mini&ci.0[i]=="black"){mini=emp.edges[i]
    }else if(emp.edges[i]>maxi&ci.0[i]=="black"){maxi=emp.edges[i]
    }else{}
  }
  
  ci.0_full=c("black","red")[((D$lower.cor>0&D$upper.cor>0)|(D$lower.cor<0&D$upper.cor<0))+1]  
  robust.indx<-which(ci.0_full=="red")
  nonrobust.indx<-which(ci.0_full=="black")
  
  return(list(thresh.min=mini,thresh.max=maxi,robust.indx=robust.indx,nonrobust.indx=nonrobust.indx))
}# End analyze.boot function

## END FUNCTIONS


# specify node shapes
# divide into groups if comparing among modalities
phenoGroupNames=rep("trait1",6)

#use 100000 for final version
pheno.boot.cor.output<-boot.cor(phenoInput,phenoGroupNames,iterations=100000)
pheno.boot.analy<-analyze.boot(pheno.boot.cor.output) #analyze and plot bootstrap output

#Get correlation matrix for unfiltered dataset
pheno.nodes<-cor(phenoInput,method="spearman")
# Remove nonrobust edges from graph
pheno.nodes[pheno.boot.analy$nonrobust.indx]<-0
diag(pheno.nodes)<-0

pheno.traitLabels=c("day","bearing","kipps","tail","tarsus","cord")

pheno.colours=c("goldenrod3","coral3","aquamarine4","aquamarine4","coral3","aquamarine4")

colnames(pheno.nodes)
pheno.out<-pheno.nodes%>%as_tibble%>%
  mutate(colourPCs=pheno.colours,rowIDs=pheno.traitLabels)%>%
  rename_at(vars(colnames(pheno.nodes)), function(x) pheno.traitLabels)
#write.csv(pheno.out,"phenotypeNetwork.csv",row.names=F)

phenoNetwork<-pheno.out
phenoNetwork<-phenoNetwork%>%rename("wing\nchord"=cord)

pheno.nodes<-phenoNetwork%>%
  mutate(rowIDs=colnames(phenoNetwork)[1:nrow(phenoNetwork)])%>%
  tibble::column_to_rownames('rowIDs')%>%
  select(1:6)%>%
  mutate_if(is.character, as.numeric)

png(filename="../figures/phenoNetwork.png",width=600,height=400)
qgraph( pheno.nodes,color="white",layout="spring",
        labels=colnames(pheno.nodes),
        minimum=0.1, #minimum edge to plot
        vsize=12, #node size
        label.norm="0000", #normalize node labels to string length
        fade=F,shape="circle",label.scale=T,
        border.color=phenoNetwork$colourPCs,border.width=8,
        negCol="grey60",posCol="grey20",label.color=1,cut=0 )

dev.off()


# Display phenotype network
Q<-qgraph( pheno.nodes,color="white",layout="spring",
           labels=pheno.traitLabels,
           minimum=0.1, #minimum edge to plot
           vsize=12, #node size
           label.norm="0000", #normalize node labels to string length
           fade=F,shape="circle",label.scale=T,
           border.color=pheno.colours,border.width=8,
           negCol="grey60",posCol="grey20",label.color=1,cut=0 )


###migratory traits

migrateGroupNames=rep("trait1",13)

migrate.boot.cor.output<-boot.cor(migrateInput,migrateGroupNames,iterations=1000)
migrate.boot.analy<-analyze.boot(migrate.boot.cor.output) #analyze and plot bootstrap output

#Get correlation matrix for unfiltered dataset
migrate.nodes<-cor(migrateInput,method="spearman")
# Remove nonrobust edges from graph
migrate.nodes[migrate.boot.analy$nonrobust.indx]<-0
diag(migrate.nodes)<-0

names(migrateInput)
migrate.traitLabels=c("day 1","day 2","day 4","bear 1","bear 2","bear 4",
              "distal","p9","p10","kipps","tail","tarsus","cord")

# Display phenotype network
Q<-qgraph( migrate.nodes,color="white",layout="spring",
           labels=migrate.traitLabels, 
           minimum=0.1, #minimum edge to plot
           vsize=10, #node size
           label.norm="00000", #normalize node labels to string length
           fade=F,shape="circle",label.scale=T,
           border.color="grey10",border.width=6,
           negCol="grey60",posCol="grey20",label.color=1,cut=0 )


###with ancestry and heterozygosity

hiGroupNames=rep("trait1",10)

#use 100000 for final version
hi.boot.cor.output<-boot.cor(hiInput,hiGroupNames,iterations=100000)
hi.boot.analy<-analyze.boot(hi.boot.cor.output) #analyze and plot bootstrap output

#Get correlation matrix for unfiltered dataset
hi.nodes<-cor(hiInput,method="spearman")
# Remove nonrobust edges from graph
hi.nodes[hi.boot.analy$nonrobust.indx]<-0
diag(hi.nodes)<-0

names(hiInput)
hi.traitLabels=c("distal","day","bearing","p9","p10","kipps","tail","tarsus","cord","ancestry")
hi.cols=c(rep("grey10",10),rep("grey10",1))

# Display phenotype network
Q<-qgraph( hi.nodes,color="white",layout="spring",
           labels=hi.traitLabels,
           minimum=0.1, #minimum edge to plot
           vsize=10, #node size
           label.norm="00000", #normalize node labels to string length
           fade=F,shape="circle",label.scale=T,
           border.color=hi.cols,border.width=6,
           negCol="grey60",posCol="grey20",label.color=1,cut=0 )


