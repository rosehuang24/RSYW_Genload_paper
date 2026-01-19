library(ggplot2)
library(reshape2)
library(dplyr)
#library(tidyverse)
library(ggpubr)
library(rstatix)
library(multcompView)

setwd("/storage/zhenyingLab/huangruoshi/genload_53/")


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
              #"long"="darkorange3",
              #"mega"="orangered4",
              "nonROH"="grey50",
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


###============
## with variant effect
# HDR
#auto cut off
###============


d<-read.delim("datatable/ROH.indv.category.class.homcount",header=TRUE)

dcoding<-d%>%filter(class %in% c("deleterious","LoF","synonymous"))
dnoncoding<-d%>%filter(class %in% c("neutral","gerp226"))

CDS<-read.delim("datatable/ROH.indv.cate.CDS.no_rep_length",header=T)
nonCDS<-read.delim("datatable/ROH.indv.cate.nonCDS.no_rep_length",header=T)

dcoding$length<-CDS[match(paste(dcoding$indv,dcoding$category),paste(CDS$indv,CDS$category)),"CDS"]
dcoding$HDR<-dcoding$count/dcoding$length

dnoncoding$length<-nonCDS[match(paste(dnoncoding$indv,dnoncoding$category),paste(nonCDS$indv,nonCDS$category)),"nonCDS"]
dnoncoding$HDR<-dnoncoding$count/dnoncoding$length


df<-dnoncoding
df<-dcoding


df<-na.omit(df)
df<-df%>%select(c("indv","category","class","HDR"))
df<-reshape(df, idvar = c("indv","class"), timevar = "category", direction = "wide")


df$`500kb-1Mb`<-df$HDR.short/df$HDR.nonROH
df$`1Mb-2Mb`<-df$HDR.medium/df$HDR.nonROH
df$`>2Mb`<-df$HDR.long/df$HDR.nonROH

df<-df%>%select(c("indv","class","500kb-1Mb","1Mb-2Mb",">2Mb"))

df<-melt(df,id=(c("indv","class")))
colnames(df)<-c("indv","class","category","ratio")
df<-na.omit(df)

df$breed<-breedinfo[match(paste(df$indv),paste(breedinfo$V2)),"V1"]



## == 
#no sep pops

stat <- df %>% 
  group_by(category) %>%
  t_test(ratio ~ class) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()%>%
  add_xy_position(x = "class")
stat
#p=0.037

ggplot(df,aes(factor(class,class_order),ratio,color=factor(class,class_order)))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(aes(shape=breed),size=1.5,width = 0.1)+
  theme_light()+
  facet_grid(category~.,scale="free",switch="y")+
  labs(x="",y="",title="Homozygotes densitiy ratio for coding variants")+
  geom_hline(yintercept=1, linetype="dashed", color="darkgrey")+
  scale_colour_manual("Variant impact",values = col.list)+
  theme(legend.title = element_text(face = "bold",size = 10),
        strip.placement = "outside",
        strip.text.y = element_text(colour = "black",face = "bold",size=16),
        strip.background = element_blank(),
        axis.text.x=element_text(colour = "black",size=13),
        axis.text.y = element_text(colour = "black",size=13),
        axis.title.y = element_text(face="bold",colour = "black",size=18),
        plot.title = element_text(face="bold",size=20),
        panel.spacing = unit(1, "lines"))+stat_pvalue_manual(stat)






#facet pops Fig 3cd:
stat <- df %>% 
  group_by(breed,category) %>%
  t_test(ratio ~ class) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()%>%
  add_xy_position(x = "class")
stat

ggplot(df,aes(factor(class,class_order),ratio,color=factor(class,class_order)))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_dodge(width = 0.75))+
  theme_light()+
  facet_grid(category~factor(breed,level_order),scale="free",switch="x")+
  labs(x="",y="short",title="Homozygotes densitiy ratio for noncoding variants")+
  geom_hline(yintercept=1, linetype="dashed", color="darkgrey")+
  scale_colour_manual("Variant impact",values = col.list)+
  theme(legend.title = element_text(face = "bold",size = 10),
        axis.text.x=element_blank(),
        axis.text.y = element_text(colour = "black",size=13),
        axis.title.y = element_text(face="bold",colour = "black",size=18),
        plot.title = element_text(face="bold",size=20),
        panel.spacing.x=unit(0, "lines"))+
  stat_pvalue_manual(stat)


