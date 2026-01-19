library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(ggpubr)
library(rstatix)
library(multcompView)
library(cowplot)


setwd("/storage/zhenyingLab/huangruoshi/genload_53/")

level_order <- c('RJF','YVC',"WHYC",'YNLC','TLF','DULO','TBC','LX','SK','WLH')
class_order<-c("LoF","deleterious","synonymous","gerp32","gerp226","gerp2","neutral")
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
              "long"="darkorange3",
              "mega"="orangered4",
              "neutral"="skyblue4", 
              "synonymous"="skyblue1", 
              "deleterious"="orange",
              "gerp2"="hotpink2",
              "gerp226"="hotpink3",
              "gerp32"="hotpink4",
              "LoF"="orangered2",
              "non_conserved"="lightskyblue2",
              "Selected"="olivedrab4",
              "Nonselected"="orange1",
              "hom"="olivedrab4",
              "het"="orange1") 






###==================================================
# DAC_SNP in the same population, in vs out, ** total sel!! -- as long as it is detected.
###==================================================

site<-read.delim("datatable/jack.breed.indv.region.class.site", header=T)
data <- read.table("datatable/jack.breed.indv.region.class.totalfreq", header = TRUE, sep = "\t")


data$site<-site[match(paste(data$jack_indv,data$region,data$class),paste(site$jack_indv,site$region,site$class)),"site"]
data$sample<-harmonic.jack[match(paste(data$breed),paste(harmonic.jack$breed)),"sample"]

data$site<-as.numeric(data$site)
data$sample<-as.numeric(data$sample)
data$DAC_SNP<-(data$totalfreq*data$sample)/data$site


d<-data%>%filter(region==breed)%>%select(c("breed","jack_indv","class","DAC_SNP"))#%>%filter(class!="synonymous")%>%filter(class!="neutral")
colnames(d)<-c("breed","jack_indv","class","Selected")
dn<-data%>%filter(region=="nonsel")
d$Nonselected<-dn[match(paste(d$jack_indv,d$class),paste(dn$jack_indv,dn$class)),"DAC_SNP"]
d<-melt(d,id=(c("breed", "jack_indv","class")))
colnames(d)<-c("breed","jack_indv","class","region","DAC_SNP")
tmp<-data[data$class=="gerp3.2",]
summary_stats <- d %>%
  group_by(class,breed,region) %>%
  summarise(
    mean_grade = mean(DAC_SNP),
    se = sd(DAC_SNP) / sqrt(n())
  )


g1<-ggplot(summary_stats[summary_stats$class == "deleterious",], aes(factor(region, levels = region_order), y=mean_grade, fill=factor(region, levels = region_order)))+
#g1<-ggplot(summary_stats[summary_stats$class == "gerp2",], aes(factor(region, levels = region_order), y=mean_grade, fill=factor(region, levels = region_order)))+
  theme_classic()+
  geom_col() +
  geom_errorbar(aes(ymin = mean_grade - se, ymax = mean_grade + se), width = 0.2) +
  scale_fill_manual(values = col.list)+
  #labs(y="synonymous\n\nDAC_SNP",x=" ")+
  labs(y="Deleterious \n\nDAC_SNP",x=" ")+
  facet_wrap(.~factor(breed, level_order),scales="free")+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black"), 
        strip.text.y = element_text(colour = "black",face = "bold",size=18, margin=margin(0,-15,0,0)),
        strip.text.x = element_text(colour = "black",face = "bold",size=16, margin=margin(0,0,15,0)),
        axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))
g1

g2<-ggplot(summary_stats[summary_stats$class == "gerp2",], aes(factor(region, levels = region_order), y=mean_grade, fill=factor(region, levels = region_order)))+
  theme_classic()+
  geom_col() +
  geom_errorbar(aes(ymin = mean_grade - se, ymax = mean_grade + se), width = 0.2) +
  scale_fill_manual(values = col.list)+
  labs(y="GERP2 \n\nDAC_SNP",x=" ",y=" ")+
  facet_wrap(.~factor(breed, level_order),scales="free")+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text  = element_blank(),
        axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))
g2

g3<-ggplot(summary_stats[summary_stats$class == "gerp226",], aes(factor(region, levels = region_order), y=mean_grade, fill=factor(region, levels = region_order)))+
  theme_classic()+
  geom_col() +
  geom_errorbar(aes(ymin = mean_grade - se, ymax = mean_grade + se), width = 0.2) +
  scale_fill_manual(values = col.list)+
  labs(y="GERP 2.26 \n\nDAC_SNP",x=" ",y=" ")+
  facet_wrap(.~factor(breed, level_order),scales="free")+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text  = element_blank(),
        axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))

g4<-ggplot(summary_stats[summary_stats$class == "gerp32",], aes(factor(region, levels = region_order), y=mean_grade, fill=factor(region, levels = region_order)))+
  theme_classic()+
  geom_col() +
  geom_errorbar(aes(ymin = mean_grade - se, ymax = mean_grade + se), width = 0.2) +
  scale_fill_manual(values = col.list)+
  labs(y="GERP 3.2\n\nDAC_SNP",x=" ")+
  facet_wrap(.~factor(breed, level_order),scales="free")+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text  = element_blank(),
        axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))

plot_grid(g1,g3,nrow=2)




###==================================================
# DAC_SNP ratio in vs out for all populations 
#***private sel region only! -- anything overlapped with the other DC is removed
###==================================================

site<-read.delim("datatable/jack.breed.indv.region.class.site.private_sel", header=T)
data <- read.table("datatable/jack.breed.indv.region.class.totalfreq.private_sel", header = TRUE, sep = "\t")

data$site<-site[match(paste(data$jack_indv,data$region,data$class),paste(site$jack_indv,site$region,site$class)),"site"]
data$sample<-harmonic.jack[match(paste(data$breed),paste(harmonic.jack$breed)),"sample"]

data$site<-as.numeric(data$site)
data$sample<-as.numeric(data$sample)
data$DAC_SNP<-(data$totalfreq*data$sample)/data$site


d<-data%>%filter(region!="nonsel")#%>%filter(class!="synonymous")%>%filter(class!="neutral")
dn<-data%>%filter(region=="nonsel")

d$nonsel<-dn[match(paste(d$jack_indv,d$class),paste(dn$jack_indv,dn$class)),"DAC_SNP"]
d$ratio<-d$DAC_SNP/d$nonsel

summary_stats <- d %>%
  group_by(class,region,breed) %>%
  summarise(
    mean_grade = mean(ratio),
    se = sd(ratio) / sqrt(n())
  )


g1<-ggplot(summary_stats[summary_stats$class == "deleterious",], 
           aes(factor(breed, levels = level_order), y=mean_grade,
               fill=factor(breed, levels = level_order)))+
  theme_light()+
  geom_col() +
  geom_text(aes(label = round(mean_grade,digits=2)),size=3, color="white", vjust = 2)+
  geom_errorbar(aes(ymin = mean_grade - se, ymax = mean_grade + se), width = 0.2) +
  scale_fill_manual(values = col.list)+
  labs(y="Deleterious\n\nDAC_SNP ratio",x=" ")+
  geom_hline(yintercept=1,linetype='dashed')+
  facet_grid(factor(class, class_order) ~ factor(region, level_order) )+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black"), 
        strip.text.y = element_text(colour = "black",face = "bold",size=18, margin=margin(0,-15,0,0)),
        strip.text.x = element_text(colour = "black",face = "bold",size=16, margin=margin(0,0,15,0)),
        axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))

g2<-ggplot(summary_stats[summary_stats$class == "gerp2",], 
           aes(factor(breed, levels = level_order), y=mean_grade,
               fill=factor(breed, levels = level_order)))+
  theme_light()+
  geom_col() +
  geom_text(aes(label = round(mean_grade,digits=2)),size=3, color="white", vjust = 2)+
  geom_errorbar(aes(ymin = mean_grade - se, ymax = mean_grade + se), width = 0.2) +
  scale_fill_manual(values = col.list)+
  labs(y="GERP2\n\nDAC_SNP ratio",x=" ")+
  geom_hline(yintercept=1,linetype='dashed')+
  facet_grid(factor(class, class_order) ~ factor(region, level_order) )+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text  = element_blank(),
        axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))
g3<-ggplot(summary_stats[summary_stats$class == "gerp226",], 
           aes(factor(breed, levels = level_order), y=mean_grade,
               fill=factor(breed, levels = level_order)))+
  theme_light()+
  geom_col() +
  geom_text(aes(label = round(mean_grade,digits=2)),size=3, color="white", vjust = 2)+
  geom_errorbar(aes(ymin = mean_grade - se, ymax = mean_grade + se), width = 0.2) +
  scale_fill_manual(values = col.list)+
  labs(y="GERP 2.26\n\nDAC_SNP ratio",x=" ")+
  geom_hline(yintercept=1,linetype='dashed')+
  facet_grid(factor(class, class_order) ~ factor(region, level_order) )+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text  = element_blank(),
        axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))

g4<-ggplot(summary_stats[summary_stats$class == "gerp32",], 
           aes(factor(breed, levels = level_order), y=mean_grade,
               fill=factor(breed, levels = level_order)))+
  theme_light()+
  geom_col() +
  geom_text(aes(label = round(mean_grade,digits=2)),size=3, color="white", vjust = 2)+
  geom_errorbar(aes(ymin = mean_grade - se, ymax = mean_grade + se), width = 0.2) +
  scale_fill_manual(values = col.list)+
  labs(y="GERP 3.2\n\nDAC_SNP ratio",x=" ")+
  geom_hline(yintercept=1,linetype='dashed')+
  facet_grid(factor(class, class_order) ~ factor(region, level_order) )+
  theme(legend.title = element_text(face = "bold",size = 8),
        axis.title.y = element_text(size=14, face="bold", colour = "black",margin=margin(0,15,0,0)), 
        strip.text  = element_blank(),
        axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=8),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines"))


#g1
plot_grid(g1,g3, nrow=2)


