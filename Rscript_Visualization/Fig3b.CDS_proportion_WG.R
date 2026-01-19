library(ggplot2)
library(reshape2)
library(dplyr)
#library(tidyverse)
library(ggpubr)
library(rstatix)
library(multcompView)

setwd("/storage/zhenyingLab/huangruoshi/genload_53/")

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

setwd("/storage/zhenyingLab/huangruoshi/genload_53/")


breedinfo<-read.table("/storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.breed.indv.txt", header=FALSE,sep="\t")
data<-read.delim("datatable/ROH.indv.category.no_rep_length",header=T)
CDSlength<-read.delim("datatable/ROH.indv.cate.CDS.no_rep_length",header=T)

data$CDS<-CDSlength[match(paste(data$indv,data$category),paste(CDSlength$indv,CDSlength$category)),"CDS"]
data$breed<-breedinfo[match(paste(data$indv),paste(breedinfo$V2)),"V1"]


data$proportion<-data$CDS/data$length
data<-na.omit(data)

cate_names <- c(
  `short` = "0.5-1Mb",
  `medium` = "1-2Mb",
  `long` = ">2Mb")


anova <- aov(proportion ~ category, data = data)
tukey <- TukeyHSD(anova)
letters <- multcompLetters4(anova, tukey)
#to output tables
table<-as.data.frame(tukey$`category`)
table <-rownames_to_column(table, var = "comparison")
write.csv(table,file="Fig3b.tukeyPs.csv",row.names = F)

#remember to detach plyr and dplyr, then re-library dplyr
Tk <- group_by(data, category) %>%
  summarise(mean = mean(proportion), 
            quant = quantile(proportion, probs = 0.75)) %>%
  mutate(letters = letters$category$Letters[category])


summary_stats <- data %>%
  group_by(category) %>%
  summarise(
    mean_grade = mean(proportion)
  )

ggplot(data,aes(factor(category,category_order),proportion,color=factor(category, category_order)))+
  geom_boxplot(outlier.colour = NA,width=0.8)+
  geom_jitter(aes(shape=breed),size=1.5,width = 0.1)+
  scale_color_manual("Length category",values = col.list)+
  theme_bw()+
  labs(x="",y="CDS Proportion")+
  scale_x_discrete(labels = c("nonROH","0.5-1Mb","1-2Mb",">2Mb"))+
  theme(legend.title = element_text(face = "bold",size = 12),
        axis.title.y = element_text(size=14, face="bold", colour = "black"), 
        strip.placement = "outside",
        strip.text.y = element_text(colour = "black",face = "bold",size=16),
       # strip.text.x = element_text(colour = "black",face = "bold",size=16),
        axis.text.x = element_text(colour = "black",face = "bold",size=14),
        axis.text.y = element_text(colour = "black",size=12),
        strip.background = element_blank(),
        panel.spacing = unit(0.8, "lines"))+
  geom_text(data = Tk, aes(x = category, y = quant, label = letters, color="black"), size = 7,
            vjust=-1,hjust=-1.5)



