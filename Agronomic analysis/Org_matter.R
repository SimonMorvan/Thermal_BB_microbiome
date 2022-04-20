# 01_02_2022
# Simon Morvan
# 

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 
#    Analysis of the soil organic matter thickness before and after burning
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
OM <- read.csv2("Org_matter.csv")

OM$Plot <- as.factor(as.character(OM$Plot))
OM$Block <- as.factor(as.character(OM$Block))
OM$Burning_intensity <- as.factor(as.character(OM$Burning_intensity))
OM$Burning_intensity <- factor(OM$Burning_intensity , levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))


####  Stats  #### 

# No burn stats 
OM_NoBurn <- subset(OM,OM$Burning_intensity=="NoBurn")
OM_NoBurn <- subset(OM_NoBurn,select=-c(Burning_intensity))
OM_NoBurn <- melt(OM_NoBurn,id.vars = c("Plot", "Block","Replicate"))
Burning_OM_NoBurn_Hist <- hist(as.numeric(OM_NoBurn$value))
Model_OM_NoBurn<- lmer(as.numeric(value) ~ variable + (1|Block), data=OM_NoBurn,REML=TRUE)
plot(fitted(Model_OM_NoBurn),residuals(Model_OM_NoBurn))
anova(Model_OM_NoBurn)

# Low burn stats 
OM_LowBurn <- subset(OM,OM$Burning_intensity=="LowBurn")
OM_LowBurn <- subset(OM_LowBurn,select=-c(Burning_intensity))
OM_LowBurn <- melt(OM_LowBurn,id.vars = c("Plot", "Block","Replicate"))
Burning_OM_LowBurn_Hist <- hist(as.numeric(OM_LowBurn$value))
Model_OM_LowBurn<- lmer(as.numeric(value) ~ variable + (1|Block), data=OM_LowBurn,REML=TRUE)
plot(fitted(Model_OM_LowBurn),residuals(Model_OM_LowBurn))
anova(Model_OM_LowBurn)


# Med burn stats 

OM_MedBurn <- subset(OM,OM$Burning_intensity=="MedBurn")
OM_MedBurn <- subset(OM_MedBurn,select=-c(Burning_intensity))
OM_MedBurn <- melt(OM_MedBurn,id.vars = c("Plot", "Block","Replicate"))
Burning_OM_MedBurn_Hist <- hist(as.numeric(OM_MedBurn$value))
Model_OM_MedBurn<- lmer(as.numeric(value) ~ variable + (1|Block), data=OM_MedBurn,REML=TRUE)
plot(fitted(Model_OM_MedBurn),residuals(Model_OM_MedBurn))
anova(Model_OM_MedBurn)

# High burn stats 
OM_HighBurn <- subset(OM,OM$Burning_intensity=="HighBurn")
OM_HighBurn <- subset(OM_HighBurn,select=-c(Burning_intensity))
OM_HighBurn <- melt(OM_HighBurn,id.vars = c("Plot", "Block","Replicate"))
Burning_OM_HighBurn_Hist <- hist(as.numeric(OM_HighBurn$value))
Model_OM_HighBurn<- lmer(as.numeric(value) ~ variable + (1|Block), data=OM_HighBurn,REML=TRUE)
plot(fitted(Model_OM_HighBurn),residuals(Model_OM_HighBurn))
anova(Model_OM_HighBurn)

# Stat table #
stat.test <- tibble::tribble(
  ~Burning_intensity, ~group1, ~group2, ~p.adj,
  "NoBurn",  "OM_before_burning", "OM_after_burning", 0.5364,
  "LowBurn",  "OM_before_burning", "OM_after_burning",0.7262 ,
  "MedBurn",  "OM_before_burning", "OM_after_burning",0.4212 ,
  "HighBurn",  "OM_before_burning", "OM_after_burning",0.4248)

stat.test


####  FIGURE S4 A  #### 
OM_long <- melt(OM,id.vars = c("Plot", "Block","Replicate","Burning_intensity"))
OM_long$value <- as.numeric(OM_long$value )
OM_long$variable 

Burning_OM<- ggplot(OM_long, aes(y=value, x=variable, colour=Burning_intensity),linetype="dotted") + 
  geom_boxplot(aes(linetype=variable))+
  geom_jitter()+
  facet_grid(.~factor(Burning_intensity, levels=c("NoBurn", "LowBurn", "MedBurn","HighBurn")),drop = TRUE, scales="free_x")+
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
  ylab("Organic layer thickness (cm)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  scale_linetype_manual(name="Status",labels=c('Before burning',"After burning"),values=c(1,2))+
  stat_pvalue_manual(stat.test,
                     y.position=6.5,
                     label="p.adj")
  
Burning_OM

#### SUP. FIGURE 4 ####

Burning_Soil_Parameters<- ggarrange(Burning_OM,
                                 Burning_humidity,
                                 labels="AUTO",
                                 ncol=1, 
                                 nrow=2, 
                                 common.legend = TRUE, 
                                 legend="bottom")
