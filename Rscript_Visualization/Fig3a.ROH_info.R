library(ggplot2)
library(reshape2)
library(dplyr)
#library(tidyverse)
library(ggpubr)
library(rstatix)
library(multcompView)

setwd("/storage/zhenyingLab/huangruoshi/genload_53/")

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- variable[,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}


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



droh <- read.table("datatable/ROH.indv.category.no_rep_length", header = TRUE)
breedinfo<-read.table("/storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.breed.indv.txt", header=FALSE,sep="\t")
droh$breed<-breedinfo[match(paste(droh$indv),paste(breedinfo$V2)),"V1"]

droh<-droh%>%filter(category!="nonROH")

droh$Froh<-droh$length/753882053

ggplot(droh,aes(factor(breed,level_order),Froh, color=factor(category, category_order)))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_dodge(width = 0.75))+
  scale_color_manual("Length category",values = col.list)+
  theme_bw()+
  labs(x="",y="FROH",title="FROH for ROHs in different length category")+
  theme(legend.title = element_text(face = "bold",size = 10),
        axis.text.x = element_text(face = "bold", colour = "black",size=14),
        axis.text.y = element_text(colour = "black",size=14),
        axis.title.y = element_text(size=14,face="bold"),
        plot.title = element_text(face="bold",size=20))

input="medium"
anova <- aov(Froh ~ breed, data = droh[droh$category==input,])
tukey <- TukeyHSD(anova)
#to output tables
filename=paste0(input, ".Fig3a.tukeyPs.csv")
table<-as.data.frame(tukey$`breed`)
table <-rownames_to_column(table, var = "comparison")
write.csv(table,file=filename,row.names = F)


letters <- multcompLetters4(anova, tukey)
letters

#short:
#WLH  SK YVC RJF 
#"a" "b" "c" "c" 
#medium
#"a"  "b" "bc"  "c" 
#long
#"a" "b" "b" "b" 



