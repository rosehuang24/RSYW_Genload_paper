library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(lemon)

setwd("/storage/zhenyingLab/huangruoshi/genload_53/")
#level_order <- c('RJF',"YVC",'SK','WLH')
level_order <- c('RJF','YVC',"WHYC",'YNLC','TLF','DULO','TBC','LX','SK','WLH')
breed_order <- c('RJF','YVC',"WHYC",'YNLC','TLF','DULO','TBC','LX','SK','WLH')

level_order <- c('RJF',"YVC",'SK','WLH')
class_order<-c("LoF","deleterious","gerp32","gerp226","gerp2","synonymous","neutral")
var_order<-c("LoF","deleterious","gerp32","gerp226","gerp2","gerp1","synonymous","neutral")
length_order<-c("all","short","medium","long","mega","nonROH")
category_order<-c("all","nonROH","short","medium","long","mega")
col.list <- c("WLH"="red", 
              "SYW"="springgreen4", 
              "SK"="Orange",
              "DULO"="deepskyblue1", 
              "YNLC"="yellow3", 
              "YVC"="violetred4", 
              "LX"="slateblue1", 
              "RJF"="grey34",
              "short"="goldenrod1",
              "long"="darkorange3",
              "mega"="orangered4",
              "neutral"="skyblue3", 
              "synonymous"="darkgoldenrod3", 
              "deleterious"="orange",
              "gerp2"="hotpink2",
              "LoF"="orangered2",
              "non_conserved"="lightskyblue2",
              "sel"="olivedrab4",
              "nonsel"="orange1",
              "hom"="olivedrab4",
              "het"="orange1") 



data <- read.table("datatable/WG.indv.zyg.class.count", header = TRUE, sep = "\t")
breedinfo<-read.table("/storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.breed.indv.txt", header=FALSE,sep="\t")
data$breed<-breedinfo[match(paste(data$indv),paste(breedinfo$V2)),"V1"]

###====================================
### === class facets on the same plot=== ###
###====================================
d<-data
d<-reshape(d, idvar=c("breed","indv","class"), timevar = c("zyg"), direction="wide")
colnames(d)<-c("indv","class","breed","Heterozygous","Homozygous")
d$Allele<-d$Homozygous*2+d$Heterozygous
d<-melt(d,id=(c("indv", "class","breed")))
colnames(d)<-c("indv","class","breed","facet","count")
label_names <- c(
#  `Homozygous` = "Homozygous \nDerived Alleles",
#  `Heterozygous` = "Heterozygous \nDerived Alleles",
#  `Allele` = "Derived Alleles"
    `Homozygous` = "Realized Load",
    `Heterozygous` = "Masked load",
    `Allele` = "Derived Alleles",
    `LoF`="LoF" ,
    `deleterious` = "Missense \nDeleterious",
    `synonymous` = "Synonymous",
    `gerp226` = "Non-coding \nDeleterious",
    `neutral` = "Neutral"
  
  )



g1<-ggplot(d[d$class == "LoF",], aes(factor(breed,breed_order),count,color = factor(class, class_order)))+
  #geom_violin(adjust = 2)+
  geom_boxplot(outlier.colour = NA)+
  theme_classic()+
  geom_point(size=1,position = position_jitter(width = 0.1))+
  scale_colour_manual(values = col.list)+
  labs(y="LoF",x=" ")+
  facet_wrap(.~factor(facet,facet_order),scales="free",labeller = as_labeller(label_names))+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black"), 
        strip.text.y = element_text(colour = "black",face = "bold",size=18, margin=margin(0,-15,0,0)),
        strip.text.x = element_text(colour = "black",face = "bold",size=16, margin=margin(0,0,15,0)),
        axis.text.x = element_text(face = "bold", colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)),labels = scales::comma)



g2<-ggplot(d[d$class == "deleterious",], aes(factor(breed,breed_order),count,color = factor(class, class_order)))+
  #geom_violin(adjust = 1.7)+
  geom_boxplot(outlier.colour = NA)+
  theme_classic()+
  geom_point(size=1,position = position_jitter(width = 0.04))+
  scale_colour_manual(values = col.list)+
  labs(y="Deleterious",title=" ",x=" ")+
  facet_wrap(.~factor(facet,facet_order),scales="free")+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text = element_blank(),
        axis.text.x = element_text(face = "bold", colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)),labels = scales::comma)

g3<-ggplot(d[d$class == "synonymous",], aes(factor(breed,breed_order),count,color = factor(class, class_order)))+
  #geom_violin(adjust = 1.7)+
  geom_boxplot(outlier.colour = NA)+
  theme_classic()+
  geom_point(size=1,position = position_jitter(width = 0.1))+
  scale_colour_manual(values = col.list)+
  labs(y="GERP2.26",title=" ",x=" ")+
  facet_wrap(.~factor(facet,facet_order),scales="free")+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text = element_blank(),
        axis.text.x = element_text(face = "bold", colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)),labels = scales::comma)


g4<-ggplot(d[d$class == "gerp226",], aes(factor(breed,breed_order),count,color = factor(class, class_order)))+
  #geom_violin(adjust = 1.7)+
  geom_boxplot(outlier.colour = NA)+
  theme_classic()+
  geom_point(size=1,position = position_jitter(width = 0.04))+
  scale_colour_manual(values = col.list)+
  labs(y="Deleterious",title=" ",x=" ")+
  facet_wrap(.~factor(facet,facet_order),scales="free")+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text = element_blank(),
        axis.text.x = element_text(face = "bold", colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)),labels = scales::comma)

g5<-ggplot(d[d$class == "neutral",], aes(factor(breed,breed_order),count,color = factor(class, class_order)))+
  #geom_violin(adjust = 1.7)+
  geom_boxplot(outlier.colour = NA)+
  theme_classic()+
  geom_point(size=1,position = position_jitter(width = 0.04))+
  scale_colour_manual(values = col.list)+
  labs(y="Deleterious",title=" ",x=" ")+
  facet_wrap(.~factor(facet,facet_order),scales="free")+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text = element_blank(),
        axis.text.x = element_text(face = "bold", colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)),labels = scales::comma)


plot_grid(g1,g2,g3,g4,g5, nrow=5)



ddd<-d%>%filter(class!="gerp2")%>%filter(class!="gerp32")
ddd<-ddd[ddd$facet!="Allele",]
g<-ggplot(ddd, aes(factor(breed,breed_order),count,color = factor(class, class_order)))+
  #geom_violin(adjust = 2)+
  geom_boxplot(outlier.colour = NA)+
  theme_bw()+
  geom_point(size=1,position = position_jitter(width = 0.1))+
  scale_colour_manual(values = col.list)+
  labs(y=" ",x=" ")+
  ggh4x::facet_grid2(factor(facet,facet_order)~factor(class,class_order),
             scales="free", independent = "y", switch="y",
             labeller = as_labeller(label_names))+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black"), 
        #strip.text.y = element_text(colour = "black",face = "bold",size=18, margin=margin(0,-15,0,0)),
        strip.text = element_text(colour = "black",face = "bold",size=16, margin=margin(0,0,15,0)),
        axis.text.x = element_text(face = "bold", colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.5, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0.4, 0.4)),labels = scales::comma)

for (facetk2 in as.character(unique(ddd$class))) {   
  for (facetk in as.character(unique(ddd$facet))) {   
    #subdf <- na.omit(subset(d, d$breed==facetk & d$class==facetk2))
    subdf<-na.omit(subset(ddd,ddd$facet==facetk & ddd$class==facetk2, select=c(breed,count)))
    model=lm(count ~ breed, data=subdf)
    ANOVA=aov(model)
    # Tukey test to study each pair of treatment :
    TUKEY <- TukeyHSD(ANOVA)
    print("The comparison is made for: ")
    #print(facetk2)
    filename=paste0(facetk2,".",facetk,".Fig2.tukeyPs.csv")
    #print(TUKEY)
    print(TUKEY$`breed`)
    table<-as.data.frame(TUKEY$`breed`)
    table <-rownames_to_column(table, var = "comparison")
    write.csv(table,file=filename,row.names = F)
    labels <- generate_label_df(TUKEY , TUKEY$`breed`)
    names(labels) <- c('Letters', 'breed')
    yvalue <- aggregate(.~breed, data=subdf,  quantile, probs=.75)  
    final <- merge(labels, yvalue)
    final$facet <-  facetk
    final$class <-  facetk2
    g <- g + geom_text(data = final,  aes(x=factor(breed, levels=level_order), y=count, label=Letters), 
                       vjust=3, hjust=1, show.legend = FALSE, size=5,color="black")}}

g

#lof or del vs syn for hom, het and all.
dms<-d%>%filter(class %in%c("deleterious","synonymous","LoF"))
dms<-reshape(dms, idvar=c("breed","indv","facet"), timevar = c("class"), direction="wide")
dms$LoFratio<-dms$count.LoF/dms$count.synonymous
dms$Delratio<-dms$count.deleterious/dms$count.synonymous
dms<-dms%>%select(!c(count.LoF,count.deleterious,count.synonymous))


#plotting
dms<-melt(dms,id=(c("indv", "facet","breed")))
colnames(dms)<-c("indv", "facet","breed","class","ratio")

g<-ggplot(dms,aes(factor(breed,level_order),ratio))+
  geom_boxplot()+
  ggh4x::facet_grid2(factor(facet,facet_order)~class, 
                     scales="free", independent = "y", switch="y")

g  
for (facetk2 in as.character(unique(dms$class))) {   
  for (facetk in as.character(unique(dms$facet))) {   
    subdf<-na.omit(subset(dms,dms$facet==facetk & dms$class==facetk2, select=c(breed,ratio)))
    model=lm(ratio ~ breed, data=subdf)
    ANOVA=aov(model)
    # Tukey test to study each pair of treatment :
    TUKEY <- TukeyHSD(ANOVA)
    print(TUKEY)
    labels <- generate_label_df(TUKEY , TUKEY$`breed`)
    names(labels) <- c('Letters', 'breed')
    yvalue <- aggregate(.~breed, data=subdf,  quantile, probs=.75)  
    final <- merge(labels, yvalue)
    final$facet <-  facetk
    final$class <-  facetk2
    g <- g + geom_text(data = final,  aes(x=factor(breed, levels=level_order), y=ratio, label=Letters), 
                       vjust=3, hjust=1, show.legend = FALSE, size=5,color="black")}}
g



#gerp vs neu for hom, het and all.
dms<-ddd%>%filter(class %in%c("gerp226","neutral"))
dms<-reshape(dms, idvar=c("breed","indv","facet"), timevar = c("class"), direction="wide")
dms$ratio<-dms$count.gerp226/dms$count.neutral
dms<-dms%>%select(!c(count.gerp226,count.neutral))
g<-ggplot(dms,aes(factor(breed,level_order),ratio))+
  geom_boxplot()+
  facet_grid(factor(facet,facet_order)~., scales="free")
g  

for (facetk in as.character(unique(dms$facet))) {   
  subdf<-na.omit(subset(dms,dms$facet==facetk , select=c(breed,ratio)))
  model=lm(ratio ~ breed, data=subdf)
  ANOVA=aov(model)
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(ANOVA)
  print(TUKEY)
  labels <- generate_label_df(TUKEY , TUKEY$`breed`)
  names(labels) <- c('Letters', 'breed')
  yvalue <- aggregate(.~breed, data=subdf,  quantile, probs=.75)  
  final <- merge(labels, yvalue)
  final$facet <-  facetk
  g <- g + geom_text(data = final,  aes(x=factor(breed, levels=level_order), y=ratio, label=Letters), 
                     vjust=3, hjust=1, show.legend = FALSE, size=5,color="black")}
g


