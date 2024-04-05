#input:
#thrush - dataframe with S, H, phenotypes, name_in_vcf
#trait.names - list of trait names, must correspond to column names in thrush dataframe
estimate_mismatch<-function(thrush,trait.names,cutoff){

  trait.l<-thrush%>%pivot_longer(cols=all_of(trait.names),
                                 values_to="trait.val",names_to="trait")%>%
    select(ancestry,heterozygosity,name.in.vcf,trait.val,trait)
  
  keep_traits<-trait.l%>%
    mutate(pop.cat=case_when(ancestry<0.25~"coastal",ancestry>0.75~"inland",TRUE~"hybrid"))%>%
    group_by(trait,pop.cat)%>%filter(!is.na(trait.val))%>%
    summarise(counts=n())%>%
    filter(pop.cat=="coastal"&counts>2)%>%pull(trait)
  trait.l<-trait.l%>%filter(trait%in%keep_traits)
  
  diff.traits.list<-find_diff_traits(trait.l,cutoff)
  
  trait.l<-trait.l%>%
    mutate("diff.trait.binom"=if_else(trait%in%diff.traits.list[[1]],"y","n"))
  
  p1<-ggplot(trait.l,aes(x=ancestry,y=trait.val,colour=diff.trait.binom))+
    geom_point(size=1,alpha=0.5)+
    geom_smooth(method="lm")+
    scale_colour_manual(values=c("grey20","firebrick"),guide="none")+
    facet_wrap(vars(trait),scales="free")+
    ylab("trait value")+xlab("")
  print(p1)
  
  hy.df<- #get values for output
    thrush%>%
    select(ancestry,heterozygosity,release.year,release.site,t5.spring40,name.in.vcf)%>%
    mutate("barrier"=est_barrier(ancestry,heterozygosity))
  
  
  trait.comp<-as.data.frame(t(combn(diff.traits.list[[1]],2)))
  colnames(trait.comp)<-c("traitname1","traitname2")
  
  
  
  #loop through trait pairs
  for(i in 1:nrow(trait.comp)){
    
    traitname1<-trait.comp[i,1]
    traitname2<-trait.comp[i,2]
    
    #reformat parents means
    m.df<-diff.traits.list[[2]]%>%
      pivot_wider(names_from="pop.type",values_from="parent.trait.val")%>%
      filter(trait%in%c(traitname1,traitname2))%>%
      mutate("trait"=c("trait1","trait2"))
    
    
    i.df<- #get hybrid values
      thrush%>%
      select(traitname1,traitname2,name.in.vcf)%>%
      drop_na()%>%
      pivot_longer(cols=c(traitname1,traitname2),values_to="tr.val",names_to="trait")%>%
      pivot_wider(values_from=tr.val,names_from=name.in.vcf)
    
    m.i<-calc_mismatch(m.df,i.df)
    
    
    i.df<-rbind(i.df,c(paste(traitname1,traitname2,sep="_"),m.i))
    i.join<-i.df%>%
      pivot_longer(!trait,names_to="name.in.vcf",values_to="tr.val")%>%
      pivot_wider(values_from=tr.val,names_from=trait)%>%
      mutate_at(vars(-name.in.vcf),as.numeric)%>%select(c(1,4))
    hy.df<-left_join(hy.df,i.join,by="name.in.vcf")
  }
  
  
  
  #mismatch for "all"

  i.df2<- #get hybrid values
    thrush%>%
    select(all_of(diff.traits.list[[1]]),name.in.vcf)%>%
    pivot_longer(cols=all_of(diff.traits.list[[1]]),values_to="tr.val",names_to="trait")%>%
    pivot_wider(values_from=tr.val,names_from=name.in.vcf)%>%
    filter(rowSums(!is.na(.))>20)%>% #keep traits with a sample size over 20
    select_if(~!any(is.na(.)))%>% #filter out individuals that have NA's
    arrange(trait)

  m.df2<-diff.traits.list[[2]]%>%
    filter(trait%in%i.df2$trait)%>%
    pivot_wider(names_from="pop.type",values_from="parent.trait.val")%>%
    arrange(trait)

  m.i2<-calc_mismatch(m.df2,i.df2)

  print(paste("mismatch traits in mismatch_all:",paste(i.df2$trait,collapse=", ")))
  i.df2<-rbind(i.df2,c("mismatch.all",m.i2))

  i.join2<-i.df2%>%
    pivot_longer(!trait,names_to="name.in.vcf",values_to="tr.val")%>%
    pivot_wider(values_from=tr.val,names_from=trait)%>%
    mutate_at(vars(-name.in.vcf),as.numeric)%>%select(c(1,mismatch.all))
  hy.df<-left_join(hy.df,i.join2)
  
  #add raw trait values to output
  hy.df<-left_join(hy.df,thrush%>%select(name.in.vcf,all_of(diff.traits.list[[2]]$trait)))
  
  
  #require a sample size of 25 for each trait or mismatch estimate
  keep.cols<-hy.df%>%
    pivot_longer(cols=c(paste(trait.comp$traitname1,trait.comp$traitname2,sep="_"),
                        diff.traits.list[[1]]),names_to="out.var")%>%
    filter(!is.na(value))%>%group_by(out.var)%>%
    summarise(count=n())%>%
    filter(count>25)%>%pull(out.var)
  hy.df<-hy.df%>%select(ancestry,heterozygosity,release.year,release.site,t5.spring40,name.in.vcf,mismatch.all,all_of(keep.cols))
  
  
  list(hy.df,diff.traits.list[[2]])}



#calculate mismatch between a set of traits
#input:
#m.df = dataframe of parent means with three cols: trait, coastal, inland and one row per trait
#i.df = dataframe of hybrid values with n cols: trait, ind1, ind2, ..., ind and one row per trait
#m.df and i.df must have the same number of rows (one per trait) with traits in the same order
#output: vector of mismatch estimates in same order as hybrid trait columns
calc_mismatch<-function(m.df,i.df){
    m.i<-c()
    for(j in 2:ncol(i.df)){
      in.co<-m.df$inland-m.df$coastal
      hy.co<-i.df%>%pull(j)-m.df$coastal
      m.j<-norm((hy.co-in.co*c(hy.co%*%in.co)/(norm(in.co,type="2")^2)),type="2")
      m.i<-c(m.i,m.j)}
    m.i}


#get list of traits with non-overlapping 75% quantiles between subspecies
#input: dataframe with trait, trait.val, and pop_type as cols
#pop_type levels should be "coastal", "inland", and "hybrid"
#trait levels should be trait names and trait.val should be numeric trait values
#output: a list of different traits
find_diff_traits<-function(trait.l,cutoff){
  
  
  diff.df<-data.frame()
  parent.preds<-data.frame()
  parents<-data.frame(ancestry=c(0,1)) #df for predicting parental values
  
  for(i in unique(trait.l$trait)){
    tr.df<-trait.l%>%filter(trait==i&!is.na(trait.val))
    lm.df<-lm(trait.val~ancestry,data=tr.df)
    aov.df<-anova(lm.df)
    #make outcome df
    diff.df<-rbind(diff.df,c(i,lm.df$coefficients[2],aov.df$`Pr(>F)`[1]))
    parent.preds<-rbind(parent.preds, #add to parent df
                        data.frame(pop.type=c("coastal","inland"),trait=i,
                                   #predict with lm and ancestry=0&1
                                   parent.trait.val=predict(lm.df,parents)))
    # 
    # tr.coastal<-tr.df%>%filter(pop.type=="coastal")%>%pull(trait.val)
    # tr.inland<-tr.df%>%filter(pop.type=="inland")%>%pull(trait.val)
    # quant.range<-max(quantile(tr.coastal)[2:4])<min(quantile(tr.inland)[2:4])|
    #   max(quantile(tr.inland)[2:4])<min(quantile(tr.coastal)[2:4])
    # if(quant.range){diff.traits<-c(diff.traits,i)}
    }
  colnames(diff.df)<-c("trait","B1","p.val")
  diff.df[,2:3]<-sapply(diff.df[,2:3],as.numeric)
  #pull out traits diff at p<0.1
  diff.traits<-diff.df%>%filter(p.val<cutoff)%>%pull(trait)
  parent.preds<-parent.preds%>%filter(trait%in%diff.traits)
  list(diff.traits,parent.preds)}

#normalize data - from Hannah
#input x is a vector of values to be normalized - no NA's
quantNorm<-function(x){qnorm(rank(x,ties.method = "average")/(length(x)+1))}


#plot pairwise mismatch results, with 
#Input: output from estimate_mismatch()
#Output: a list of plots, with last plot being the legend
plot_pairwise_mismatch<-function(mis.thrush){
  #mis.thrush<-RC.mis
  
  parents.w<-mis.thrush[[2]]%>%pivot_wider(names_from=trait,values_from=parent.trait.val)
  mis.df<-mis.thrush[[1]]
  #colnames(mis.df)<-gsub(".corr","",gsub("qnorm.","",colnames(mis.df)))
  colnames(mis.df)<-gsub("\\+",".",gsub(" ",".",colnames(mis.df)))
  colnames(parents.w)<-gsub("\\+",".",gsub(" ",".",colnames(parents.w)))
  mis.plots<-list()
  cols.i<-grep("_",colnames(mis.df)) 
  mis.df$t5.spring40<-as.character(mis.df$t5.spring40)
  
  for(i in 1:length(cols.i)){
    
    
    #scaling for viridis
    vir.max.i<-(max(mis.df[,cols.i[i]],na.rm=T)-min(mis.df[,cols.i],na.rm=T))/
      (max(mis.df[,cols.i],na.rm=T)-min(mis.df[,cols.i],na.rm=T))
    vir.min.i<-(min(mis.df[,cols.i[i]],na.rm=T)-min(mis.df[,cols.i],na.rm=T))/
      (max(mis.df[,cols.i],na.rm=T)-min(mis.df[,cols.i],na.rm=T))
    
    mis.pair<-colnames(mis.df)[cols.i][i]
    
    x.lab.1<-gsub(".corr","",gsub("qnorm.","",str_split(mis.pair,"_")[[1]][1]))%>%
      dplyr::recode("fall.bear.1"="fall bearing","doy.fall.r1"="fall departure day",
             "tarsus.length"="tarsus (leg)","PS.kipps"="PS kipps (wing)",
             "p9"="P9 (wing)","wing.cord"="wing chord")
    y.lab.1<-gsub(".corr","",gsub("qnorm.","",str_split(mis.pair,"_")[[1]][2]))%>%
      dplyr::recode("fall.bear.1"="fall bearing","doy.fall.r1"="fall departure day",
             "tarsus.length"="tarsus (leg)","PS.kipps"="PS kipps (wing)",
             "p9"="P9 (wing)","wing.cord"="wing chord")
    
    
    #pal<-rev(natparks.pals("CapitolReef",n=100,type="continuous")[1:80])
    
    mis.plots[[i]]<-
      ggplot()+
      geom_point(data=mis.df,
                 aes_string(x=str_split(mis.pair,"_")[[1]][1],y=str_split(mis.pair,"_")[[1]][2],
                            colour=mis.pair),size=2)+
      geom_point(data=parents.w,aes_string(x=str_split(mis.pair,"_")[[1]][1],
                                           y=str_split(mis.pair,"_")[[1]][2],
                                           shape="pop.type"),size=5,stroke=1.5,
                 colour="grey60")+
      scale_shape_manual(values=c(2,5),name="")+
      scale_alpha(range=c(0.3,1))+
      scale_colour_viridis(name="mismatch",option="F",begin=vir.min.i,end=vir.max.i)+
      #scale_colour_gradientn(colours=pal)+
      xlab(gsub("\\."," ",gsub("\\.\\.\\."," + ",x.lab.1)))+
      ylab(gsub("\\."," ",gsub("\\.\\.\\."," + ",y.lab.1)))+
      theme(legend.position="none",aspect.ratio=1)
    
    }
  
  #long data to get full range of values for legend
  mis.legend<-mis.df%>%pivot_longer(cols=cols.i,names_to="mismatch.pairs",values_to="mismatch.vals")
  #pal<-rev(natparks.pals("CapitolReef",n=100,type="continuous")[1:80])
  #pal<-wesanderson::wes_palette("Royal2",n=100,type="continuous")[1:35]
  l1<-ggplot()+
    geom_point(data=mis.legend,aes_string(x="ancestry",y="heterozygosity",colour="mismatch.vals"))+
    geom_point(data=parents.w,aes_string(x=str_split(mis.pair,"_")[[1]][1],
                                         y=str_split(mis.pair,"_")[[1]][2],
                                         shape="pop.type"),size=5,stroke=1.5)+
    scale_shape_manual(values=c(2,5),name="")+
    scale_colour_viridis(name="mismatch",option="rocket")+
    #scale_colour_gradientn(colours=pal,name="mismatch")+
    theme(legend.box="horizontal")
  
  mis.plots[[i+1]]<-cowplot::get_legend(l1)
  
  mis.plots}


#convert predictSurface() object to something plottable in ggplot
#input:
# HI_df=dataframe to feed into Tps, with columns ancestry, heterozygosity, and retrieved_spring
make_plot_surface<-function(HI_df,survival_var,ancestry_var,heterozygosity_var){
  # HI_df<-HI_juvie
  # survival_var<-"t5_spring40"
  # ancestry_var="ancestry"
  # heterozygosity_var="heterozygosity"
  # 
  
  tps.pop<-Tps(as.matrix(HI_df%>%select(all_of(ancestry_var),all_of(heterozygosity_var))),
               HI_df%>%select(all_of(survival_var)),method="REML")
  #surface(tps.pop)
  pred.pop<-predictSurface(tps.pop,nx=100,ny=100)
  
  z.sample<-as_tibble(data.frame(pred.pop$z))%>%
    gather()%>%
    #use existing Y values, adjust to fall between observed values
    mutate(Yval=sort(rep(pred.pop$y,100)))%>%
    #get sequence of x values that fall between observed values
    mutate(Xval=rep(pred.pop$x,100))%>%
    dplyr::rename(Zval=value)%>%
    select(Xval,Yval,Zval)%>%
    drop_na
  
  
  z.sample}

#input: dataframe that includes columns named "ancestry" and "heterozygosity"
#output: dataframe of coordinates for region of the triangle to cover in figure 
##with a polygon because there is no data there
whiteOut_tri<-function(HI_df,ancestry_var,heterozygosity_var){
  # HI_df=HI_juvie
  # ancestry_var="ancestry_scaf4"
  # heterozygosity_var="heterozygosity_scaf4"
  
  #convert character variables to symbol for filtering
  #then use !! before variable names to inject them into dplyr function
  heterozygosity_var1=as.symbol(heterozygosity_var)
  ancestry_var1=as.symbol(ancestry_var)
  
  #get max ancestry on the coastal side (bottom left of polygon)
  min_ancestry<-max(HI_df%>%filter(!!heterozygosity_var1<0.25&!!ancestry_var1<0.5)%>%
                      pull(!!ancestry_var1))+0.05
  #min ancestry on inland side (bottom right)
  max_ancestry<-min(HI_df%>%filter(!!heterozygosity_var1<0.25&!!ancestry_var1>0.5)%>%
                      pull(!!ancestry_var1))-0.05
  #min heterozygosity for F1ish hybrids - top y coord, with 0.47 and 0.53 as x coords
  max_heterozygosity<-min(HI_df%>%filter(!!ancestry_var1>0.4&!!ancestry_var1<0.6)%>%
                            pull(!!heterozygosity_var1))-0.05
  
  poly.out<-data.frame(X.co=c(min_ancestry,0.47,0.53,max_ancestry),
                       Y.co=rep(c(0,max_heterozygosity,max_heterozygosity,0)))
  
  poly.out}

