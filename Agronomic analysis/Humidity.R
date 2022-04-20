# 01_02_2022
# Simon Morvan
# 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 
#          Analysis of the soil humidity before and after burning
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#



#### PACKAGES ####
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lmerTest)
library(reshape2)

#### DATA ####

setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/Data/Soil/")
Humid <- read.csv2("Humidity.csv")

Humid$Plot <- as.factor(as.character(Humid$Plot))
Humid$Block <- as.factor(as.character(Humid$Block))
Humid$Burning_intensity <- as.factor(as.character(Humid$Burning_intensity))
Humid$Burning_intensity <- factor(Humid$Burning_intensity , levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))


####  Stats  #### 

# No burn stats 
Humid_NoBurn <- subset(Humid,Humid$Burning_intensity=="NoBurn")
Humid_NoBurn <- subset(Humid_NoBurn,select=-c(Burning_intensity))
Humid_NoBurn <- melt(Humid_NoBurn,id.vars = c("Plot", "Block","Replicate"))
Burning_humidity_NoBurn_Hist <- hist(as.numeric(Humid_NoBurn$value))
Model_humidity_NoBurn<- lmer(as.numeric(value) ~ variable + (1|Block), data=Humid_NoBurn,REML=TRUE)
plot(fitted(Model_humidity_NoBurn),residuals(Model_humidity_NoBurn))
anova(Model_humidity_NoBurn)

# Low burn stats 
Humid_LowBurn <- subset(Humid,Humid$Burning_intensity=="LowBurn")
Humid_LowBurn <- subset(Humid_LowBurn,select=-c(Burning_intensity))
Humid_LowBurn <- melt(Humid_LowBurn,id.vars = c("Plot", "Block","Replicate"))
Burning_humidity_LowBurn_Hist <- hist(as.numeric(Humid_LowBurn$value))
Model_humidity_LowBurn<- lmer(as.numeric(value) ~ variable + (1|Block), data=Humid_LowBurn,REML=TRUE)
plot(fitted(Model_humidity_LowBurn),residuals(Model_humidity_LowBurn))
anova(Model_humidity_LowBurn)


# Med burn stats 

Humid_MedBurn <- subset(Humid,Humid$Burning_intensity=="MedBurn")
Humid_MedBurn <- subset(Humid_MedBurn,select=-c(Burning_intensity))
Humid_MedBurn <- melt(Humid_MedBurn,id.vars = c("Plot", "Block","Replicate"))
Burning_humidity_MedBurn_Hist <- hist(as.numeric(Humid_MedBurn$value))
Model_humidity_MedBurn<- lmer(as.numeric(value) ~ variable + (1|Block), data=Humid_MedBurn,REML=TRUE)
plot(fitted(Model_humidity_MedBurn),residuals(Model_humidity_MedBurn))
anova(Model_humidity_MedBurn)

# High burn stats 
Humid_HighBurn <- subset(Humid,Humid$Burning_intensity=="HighBurn")
Humid_HighBurn <- subset(Humid_HighBurn,select=-c(Burning_intensity))
Humid_HighBurn <- melt(Humid_HighBurn,id.vars = c("Plot", "Block","Replicate"))
Burning_humidity_HighBurn_Hist <- hist(as.numeric(Humid_HighBurn$value))
Model_humidity_HighBurn<- lmer(as.numeric(value) ~ variable + (1|Block), data=Humid_HighBurn,REML=TRUE)
plot(fitted(Model_humidity_HighBurn),residuals(Model_humidity_HighBurn))
anova(Model_humidity_HighBurn)

# Stat table #
stat.test <- tibble::tribble(
  ~Burning_intensity, ~group1, ~group2,   ~p.adj,
     "NoBurn",  "Humidity_before_burn", "Humidity_after_burn", 0.9407,
     "LowBurn",  "Humidity_before_burn", "Humidity_after_burn",0.6683 ,
     "MedBurn",  "Humidity_before_burn", "Humidity_after_burn",0.4078 ,
     "HighBurn",  "Humidity_before_burn", "Humidity_after_burn",0.601)
  
stat.test


####  FIGURE S4 B #### 
humid_long <- melt(Humid,id.vars = c("Plot", "Block","Replicate","Burning_intensity"))
humid_long$value <- as.numeric(humid_long$value )

Burning_humidity <- ggplot(humid_long, aes(y=value, x=variable, colour=Burning_intensity),linetype="dotted") + 
  geom_boxplot(aes(linetype=variable))+
  geom_jitter()+
  facet_grid(.~factor(Burning_intensity, levels=c("NoBurn", "LowBurn", "MedBurn","HighBurn")),drop = T,scales="free_x")+
  stat_summary(fun="mean", geom="point",size=3,aes(x=variable,y=value),
               position=position_dodge(width=0.75),colour="black")+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=2,aes(x=variable,y=value,colour=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"),
        legend.box="vertical")+
  ylab(expression(paste("Humidity content ", (m^3/m^3))))+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  scale_linetype_manual(name="Status",labels=c('Before burning',"After burning"),values=c(1,2))+
  stat_pvalue_manual(stat.test,
                     y.position=0.33,
                     label="p.adj")
