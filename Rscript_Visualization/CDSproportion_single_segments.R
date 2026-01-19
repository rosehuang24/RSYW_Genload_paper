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


CDS_order<-c("all","nonCDS")
#process input files
data <- read.delim("mid_process_file/ROH/single_segments/indv.segment.repROH.norepROH.CDS", header = TRUE, sep = "\t",row.names=NULL)

breedinfo<-read.table("/storage/zhenyingLab/huangruoshi/txt_might_be_useful/RSYW.breed.indv.txt", header=FALSE,sep="\t")
data$breed<-breedinfo[match(paste(data$indv),paste(breedinfo$V2)),"V1"]
data$proportion<-data$CDS/data$norepROH
data$rep<-data$repROH-data$norepROH

bin_size <- 200000
min_score <- 100000
#max_score <- 7400000
max_score <- 5000000
bins <- c(seq(min_score, max_score, by = bin_size), Inf)
bins <- c(100000, 300000,500000, 1000000,2000000,Inf)

dd <- data %>%
  mutate(bin = cut(repROH, breaks = bins, labels = bins[-length(bins)], include.lowest = TRUE))
  
dd<-dd%>%
  group_by(indv,bin) %>% 
  summarise(across(c(norepROH,CDS), sum))

dd$proportion<-dd$CDS/dd$norepROH
dd$breed<-breedinfo[match(paste(dd$indv),paste(breedinfo$V2)),"V1"]

ggplot(dd,aes(bin,proportion))+
  geom_boxplot(outlier.colour = NA,width=0.65)+
  geom_jitter(aes(color=breed),size=1,width = 0.08)+
  scale_colour_manual(values = col.list)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_hline(yintercept=0.034, linetype="dashed", color="darkgrey")

anova <- aov(proportion ~ bin, data = dd)
tukey <- TukeyHSD(anova)
letters <- multcompLetters4(anova, tukey)
Tk <- group_by(dd, bin) %>%
  summarise(mean = mean(proportion), 
            quant = quantile(proportion, probs = 0.75)) %>%
  mutate(letters = letters$bin$Letters[bin])
Tk

#95CI
ddd <- dd%>%
  group_by(bin) %>%
  dplyr::summarise(mean = mean(proportion), 
                   upper = CI(proportion)[1], 
                   lower = CI(proportion)[3])

ggplot(ddd, aes(x=bin, y=mean))+
  geom_point(size=3, position=position_dodge(width=0.4))+
  theme_bw()+
  geom_hline(yintercept=0.034, linetype="dashed", color="darkgrey")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2,position=position_dodge(width=0.4))



d<-data%>%filter(data$norepROH>500000)
ggplot(d[d$breed=="WLH",], aes(norepROH,proportion))+
  geom_point(alpha=0.3)+
  scale_x_log10() +
  #scale_colour_manual(values = col.list)+
  geom_smooth(method = "lm")


d.lm = lm(proportion ~ repROH, data=data)
summary(d.lm)



ggplot(data[data$ROH!=0,], aes(ROH,CDS))+
  geom_point(alpha=0.5)+
  #scale_colour_manual(values = col.list)+
  geom_smooth(method = "lm")+
  geom_abline(intercept = 0, slope = 0.0347, color="red", 
              linetype="dashed", size=1.5)


d<-data[data$ROH!=0,]
model <- lm(rep ~ repROH, data = data)

summary(model)

slope <- coef(model)["ROH"]
se <- summary(model)$coefficients["ROH", "Std. Error"]
t_stat <- (slope - 0.0347) / se
d_res <- df.residual(model)
p_val <- 2 * pt(-abs(t_stat), d_res)

p_val

#distribution of ROH length#

bin_size <- 100000
min_score <- 0
#max_score <- 7400000
max_score <- 3200000
bins <- seq(min_score, max_score, by = bin_size)

data<-read.delim("ROH_sheets/input4.RSYW50.ROHs", col.names = c("pop","indv","chrom", "start", "end"))
data$plink<-data$end-data$start
garlic<-read.delim("GARLIC_ROH/man_cutoff/indv.category.no_rep.length")
dd <- data %>%
  mutate(bin = cut(ROH, breaks = bins, labels = bins[-length(bins)], include.lowest = TRUE))
counts <- dd %>%
  group_by(bin) %>%
  summarise(count = n())

ggplot(counts, aes(x = as.numeric(as.character(bin)), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  #labs(x = "GERP Score", y = "Count") +
  #ggtitle("GERP Score Distribution") +
  scale_y_log10() +
  theme_minimal()+
  scale_x_continuous(breaks = seq(min_score, max_score, by = 100000))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




data<-read.delim("../chicken_ref/repetitive_auto.bed", col.names = c("chrm", "start", "end",))
data$rep_length<-data$end-data$start
ggplot(data,aes(as.numeric(rep_length)))+
  geom_histogram(bins = 100)+
  scale_x_log10() 


dd <- data %>%
  mutate(bin = cut(ROH, breaks = bins, labels = bins[-length(bins)], include.lowest = TRUE))
counts <- dd %>%
  group_by(bin) %>%
  summarise(count = n())


ggplot(counts, aes(x = as.numeric(as.character(bin)), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  #labs(x = "GERP Score", y = "Count") +
  #ggtitle("GERP Score Distribution") +
  scale_y_log10() +
  theme_minimal()+
  scale_x_continuous(breaks = seq(min_score, max_score, by = 100000))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


