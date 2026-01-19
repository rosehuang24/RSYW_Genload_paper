library(ggplot2)
library(reshape2)
library(dplyr)
#library(tidyverse)
library(ggpubr)
library(rstatix)
library(multcompView)

setwd("/storage/zhenyingLab/huangruoshi/genload/")

#rep-free autosomal length: 753882071 ($REF/auto.no_rep.bed)
#rep-free CDS length: 25342528 (neutral_region-work/gff3.release110.CDS.no_rep.bed)
#proporion of CDS in WG: 0.03361604

level_order <- c('RJF',"YVC",'SK','WLH')
class_order<-c("LoF","deleterious","gerp32","gerp226","gerp2","synonymous","neutral")
var_order<-c("LoF","deleterious","gerp32","gerp226","gerp2","gerp1","synonymous","neutral")
length_order<-c("all","short","medium","long","mega","nonROH")
category_order<-c("all","nonROH","short","medium","long","mega")
#category_order<-c("short(100kb~300kb)","long(300kb~1Mb)","mega(>1Mb)","nonROH")
col.list <- c("WLH"="red", 
              "SYW"="springgreen4", 
              "SK"="Orange",
              "DULO"="deepskyblue1", 
              "YNLC"="yellow3", 
              "YVC"="violetred4", 
              "LX"="slateblue1", 
              "RJF"="grey34",
              "short"="goldenrod1",
              "medium"="darkorange3",
              "long"="orangered4",
              "neutral"="skyblue4", 
              "synonymous"="skyblue1", 
              "deleterious"="orange",
              "gerp2"="hotpink2",
              "gerp226"="hotpink3",
              "gerp32"="hotpink4",
              "LoF"="orangered2",
              "non_conserved"="lightskyblue2",
              "sel"="olivedrab4",
              "nonsel"="orange1",
              "hom"="olivedrab4",
              "het"="orange1") 

#setwd("/storage/zhenyingLab/huangruoshi/genload_53/ROH_sheets/chrm1-5")


breedinfo<-read.table("/storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.breed.indv.txt", header=FALSE,sep="\t")
data<-read.delim("datatable/ROH.indv.category.cate.chrm1-5.norep.length",header=T)
CDSlength<-read.delim("datatable/ROH.indv.category.CDS.chrm1-5.norep.length",header=T)

data$CDS<-CDSlength[match(paste(data$indv,data$category,"chrm",data$chrm),paste(CDSlength$indv,CDSlength$category,"chrm",CDSlength$chrm)),"CDS"]
data$breed<-breedinfo[match(paste(data$indv),paste(breedinfo$V2)),"V1"]


data$proportion<-data$CDS/data$length
data<-na.omit(data)

label_names <- c(
  `1` = "Chr 1",
  `2` = "Chr 2",
  `3` = "Chr 3",
  `4` = "Chr 4",
  `5` = "Chr 5")

cate_names <- c(
  `short` = "500kb-1Mb",
  `medium` = "1Mb-2Mb",
  `long` = ">2Mb")

#[!data$chrm %in% c("1", '2',"3"),]

g<-ggplot(data,aes(factor(category,category_order),proportion,color=factor(category, category_order)))+
  geom_boxplot(outlier.colour = NA,width=0.8)+
  geom_jitter(aes(shape=breed),size=1.5,width = 0.1)+
  scale_color_manual("Length category",values = col.list)+
  theme_bw()+
  labs(x="",y="")+
  # stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) + 
  facet_grid(chrm~., scales="free",labeller = as_labeller(label_names), switch = "y")+
  #facet_grid(chrm~., scales="free", switch = "y")+
  scale_x_discrete(labels = c("nonROH","0.5-1Mb","1Mb-2Mb",">2Mb"))+
  theme(legend.title = element_text(face = "bold",size = 12),
        axis.title.y = element_text(size=14, face="bold", colour = "black"), 
        strip.placement = "outside",
        strip.text.y = element_text(colour = "black",face = "bold",size=16),
       # strip.text.x = element_text(colour = "black",face = "bold",size=16),
        axis.text.x = element_text(colour = "black",face = "bold",size=14),
        axis.text.y = element_text(colour = "black",size=12),
        strip.background = element_blank(),
        panel.spacing = unit(0.8, "lines"))

g

for (facetk in as.character(unique(data$chrm))) {   
  subdf<-na.omit(subset(data,data$chrm==facetk, select=c(category,proportion)))
  model=lm(proportion ~ category, data=subdf)
  ANOVA=aov(model)
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(ANOVA)
  filename=paste0("chrm.",facetk,".FigS2.tukeyPs.csv")
  #print(TUKEY)
  print(TUKEY$`category`)
  table<-as.data.frame(TUKEY$`category`)
  table <-rownames_to_column(table, var = "comparison")
  write.csv(table,file=filename,row.names = F)
  labels <- generate_label_df(TUKEY , TUKEY$`category`)
  names(labels) <- c('Letters', 'category')
  yvalue <- aggregate(.~category, data=subdf,  quantile, probs=.5)  
  final <- merge(labels, yvalue)
  final$chrm <-  facetk
  print(final)
  g <- g + geom_text(data = final,  aes(x=factor(category, levels=category_order), y=proportion, label=Letters), 
                     vjust=-2, hjust=-2, show.legend = FALSE, size=4,color="black")
}

g




#percentage of chrm1-5 ROHs to WG
data<-read.delim("datatable/ROH.indv.category.cate.chrm1-5.norep.length",header=T)
droh <- read.table("datatable/ROH.indv.category.no_rep_length", header = TRUE)

sdata <- group_by(data,indv,category) %>%
  summarise(sum5chrm = sum(length))

sdata$totalROH<-droh[match(paste(sdata$indv,sdata$category),paste(droh$indv,droh$category)),"length"]
sdata$proportion<-sdata$sum5chrm/sdata$totalROH
sdata<-na.omit(sdata)

pdata <- group_by(sdata,category) %>%
  summarise(mean = mean(proportion)*100,
            se = 100*sd(proportion) / sqrt(n()))

