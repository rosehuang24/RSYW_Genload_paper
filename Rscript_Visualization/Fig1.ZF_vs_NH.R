library(ggplot2)
library(reshape2)
library(dplyr)
#library(tidyverse)
library(ggpubr)
library(rstatix)
library(multcompView)

setwd("/storage/zhenyingLab/huangruoshi/genload/")
#level_order <- c('RJF',"YVC",'SK','WLH')
level_order <- c('RJF','YVC',"WHYC",'YNLC','TLF','DL','TBC','LX','SK','WLH')
class_order<-c("high","deleterious","gerp2","gerp1","synonymous","neutral")
#var_order<-c("high","deleterious","gerp2","gerp1","synonymous","neutral")
length_order<-c("short","long","mega","nonROH")
category_order<-c("short(100kb~300kb)","long(300kb~1Mb)","mega(>1Mb)","nonROH")
col.list <- c("WLH"="red", 
              "TLF"="plum2", 
              "LX"="tan",
              "WHYC"= "palegreen1", 
              "DL"="darkseagreen", 
              "YNLC"="thistle1", 
              "SK"="lightblue", 
              "YVC"="khaki", 
              "TBC"="gray88", 
              "RJF"="blue",
              "short"="goldenrod1",
              "long"="darkorange3",
              "mega"="orangered4",
              "neutral"="skyblue3", 
              "synonymous"="darkgoldenrod3", 
              "deleterious"="orangered1",
              "gerp2"="tomato2",
              "gerp1"="tomato4",
              "high"="pink",
              "non_conserved"="lightskyblue2",
              "sel"="olivedrab4",
              "nonsel"="orange1",
              "hom"="olivedrab4",
              "het"="orange1") 
#awk '{sum+=$3-$2} END {print sum}' neutral_region-work/non_CDS_flanking_no_rep.bed
#174348285

data <- read.table("/storage/zhenyingLab/huangruoshi/108/zerofour/dataframe", header = TRUE, sep = "\t")
breedinfo<-read.table("/storage/zhenyingLab/huangruoshi/txt_might_be_useful/107.breed_indv_depth.DL.txt", header=FALSE,sep="\t")
data$breed<-breedinfo[match(paste(data$indv),paste(breedinfo$V2)),"V1"]
#data$size <- ifelse(data$breed %in% c("WLH", "RJF"), 1.8, 1.79)
data$NH<-data$neutral/174348285


ggplot(data,aes(NH,ZF))+
  geom_point(aes(colour=factor(breed,level=level_order)))+
  scale_colour_manual("Breed",values = col.list)+
  geom_smooth(method = 'lm', se=FALSE, formula = y ~ x, color="grey45",linewidth=0.7)+
  theme_classic() +
  theme(legend.title = element_text(face = "bold",size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=16, face="bold")) +
  scale_x_continuous(labels = scales::comma)+
  labs(y="0/4-fold Heterozygosity", x= "Neutral Heterozygosity")

d.lm = lm(ZF ~ NH, data=data)
summary(d.lm)

