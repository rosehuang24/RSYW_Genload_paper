library(ggplot2)
library(reshape2)
library(dplyr)
#library(tidyverse)
library(ggpubr)
library(rstatix)
library(multcompView)

setwd("/storage/zhenyingLab/huangruoshi/genload_53/")


#####===========================================
# Jckknifing blocks for wholge genome, all pairs
#####===========================================
jackknife_ratio <- function(rx, ry, k = 100) {
  n <- length(rx)
  indices <- sample(rep(1:k, length.out = n))
  jack_vals <- numeric(k)
  
  for (i in 1:k) {
    rx_sum <- sum(rx[indices != i])
    ry_sum <- sum(ry[indices != i])
    jack_vals[i] <- rx_sum / ry_sum
  }
  
  jack_mean <- mean(jack_vals)
  jack_se <- sqrt((k - 1) / k * sum((jack_vals - jack_mean)^2))
  
  return(list(mean = jack_mean, se = jack_se))
}

# List all the 20 files in your working directory
file_list <- list.files(pattern = "*.bed")


# Initialize result dataframe
results <- data.frame()


# Define populations
pops <- c("RJF", "YVC", "SK", "WLH")

setwd("/storage/zhenyingLab/huangruoshi/genload_53/datatable/RXY/")
# Iterate through each file
for (file in file_list) {
  
  # Extract region and class from filename
  match <- str_match(file, "ALT\\.freq\\.RSYW\\.(.*?)\\.bed")

  class <- match[2]
  
  # Skip unmatched files
  if (is.na(class)) next
  
  # Read the data
  df <- read.table(file, header = FALSE)
  colnames(df) <- c("A", "B", "C", "RJF", "YVC", "SK", "WLH")
  
  # Compute Rxy for every ordered pair of distinct populations
  for (popx in pops) {
    for (popy in pops) {
      if (popx == popy) next
      
      rx <- df[[popx]] * (1 - df[[popy]])
      ry <- df[[popy]] * (1 - df[[popx]])
      
      # Filter for valid data
      valid <- is.finite(rx) & is.finite(ry)
      rx <- rx[valid]
      ry <- ry[valid]
      
      if (length(rx) == 0) next
      
      jack <- jackknife_ratio(rx, ry)
      
      results <- rbind(results, data.frame(
        #file_region = region_raw,
        class = class,
        Rx = popx,
        Ry = popy,
        #rx_sum = sum(rx),
       # ry_sum = sum(ry),
        jack_mean = jack$mean,
        jack_se = jack$se
      ))
    }
  }
}

results$pairs<-paste(results$Rx,results$Ry, sep = "_")

results$rxy_shifted <- results$jack_mean - 1


results$lower<-results$rxy_shifted-results$jack_se#*1.96
results$upper<-results$rxy_shifted+results$jack_se#*1.96
results<-results%>%filter(class!="gerp1")%>%filter(class!="high")
 
pair_order<-c("RJF_YVC","RJF_SK","RJF_WLH","YVC_SK","YVC_WLH")
pair_order<-c("YVC_RJF","SK_RJF","WLH_RJF")
r<-results%>%filter(pairs %in% pair_order)%>%filter(class!="gerp2" & class!="gerp32")


coding<-c("LoF","deleterious","synonymous")
noncoding<-c("neutral","gerp226")
r<-r%>%mutate(code = case_when(class %in% coding ~ "coding",
                           class %in% noncoding ~"noncoding"))

ggplot(r, aes(x=factor(pairs,pair_order), y=rxy_shifted,width=0.8, fill=factor(class, levels=rev(class_order))))+
  #geom_point(size=3, position=position_dodge(width=0.4))+
  geom_bar(stat = "identity", position = "dodge",width=0.8)+
  theme_light()+
  facet_grid(factor(pairs,pair_order)~code, scales="free")+
  scale_fill_manual("Variant impact",values = col.list)+
  labs(x="",y="",title="Rxy")+
  theme(legend.title = element_text(face = "bold",size = 10),
        axis.text=element_text(size=12,face="bold"),
        panel.spacing.y = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face="bold",size=20))+
  scale_y_continuous(labels = function(x) x + 1,name = "Rxy")+
  scale_x_discrete(limits=rev)+
  geom_hline(yintercept=0, linetype="dashed", color="darkgrey")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2,position=position_dodge(width=.8))+
  coord_flip()

#no neutral or syn
r2<-r[r$class!="synonymous" & r$class!="neutral",]%>%filter(!pairs %in% c("YVC_SK","YVC_WLH"))
ggplot(r2, aes(x=factor(pairs,pair_order), y=rxy_shifted,width=0.8, fill=factor(class, levels=rev(class_order))))+
  #geom_point(size=3, position=position_dodge(width=0.4))+
  geom_bar(stat = "identity", position = "dodge",width=2)+
  theme_light()+
  facet_grid(factor(pairs,pair_order)~., scales="free")+
  scale_fill_manual("Variant impact",values = col.list)+
  labs(x="",y="",title="Rxy")+
  theme(legend.title = element_text(face = "bold",size = 10),
        axis.text=element_text(size=12,face="bold"),
        panel.spacing.y = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face="bold",size=20))+
  scale_y_continuous(labels = function(x) x + 1,name = "Rxy")+
  scale_x_discrete(limits=rev)+
  geom_hline(yintercept=0, linetype="dashed", color="darkgrey")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2,position=position_dodge(width=.8))+
  coord_flip()


