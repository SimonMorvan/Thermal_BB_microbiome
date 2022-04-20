# 19_10_2021
# Simon Morvan
# 

R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 
# Analysis of the blueberry coverage and biomass as well as soil coverage
#
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# Packages 
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(lme4)
library(lmerTest)
library(multcomp)

# Data 
setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/Data/Bleuet/")
BB_data<- read.csv2("BB_coverage_biomass.csv")
BB_data$Burning_intensity<- factor(BB_data$Burning_intensity , levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))
BB_data$Date <- as.factor(BB_data$Date)
BB_data$Date<- factor(BB_data$Date, levels = c("2018/06/07","2018/06/26","2018/07/17","2018/08/16","2019/09/06"))


####______________####
###### RAW DATA DISTRIBUTION #####

####______________####
#### Soil coverage #### 
BB_data$Soil_coverage <- round(BB_data$Soil_coverage,digits = 0)
BB_data$Soil_coverage <- as.integer(BB_data$Soil_coverage)
Soil_hist <- hist(BB_data$Soil_coverage)

#### Blueberry coverage #### 
BB_data$BB_coverage<- round(BB_data$BB_coverage,digits = 0)
BB_data$BB_coverage<- as.integer(BB_data$BB_coverage)
Blueberry_hist <- hist(BB_data$BB_coverage)

#### Blueberry biomass  #### 

BlueberryBiomass_hist <- hist(BB_data$BB_biomass)


####______________####
###### MODELS & BOXPLOTS #####

####______________####

####  Soil coverage #####


BB_data$Soil_coverage_TR <- asin(sqrt((BB_data$Soil_coverage/100)))
hist(BB_data$Soil_coverage_TR)

BB_data$Soil_coverage_sqrt <- sqrt(BB_data$Soil_coverage)
hist(BB_data$Soil_coverage_sqrt)


model_SoilCover = lmer(Soil_coverage_sqrt ~ Burning_intensity * Date + (1|Block/Plot), data=BB_data,REML=TRUE)

plot(fitted(model_SoilCover),residuals(model_SoilCover)) # Needs to have no particular pattern

# I. Linearity assumption 
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# III. Homoskedasticity
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_SoilCover))
# Histogram has to be bell shaped
qqnorm(residuals(model_SoilCover))
# Straight line = normal distribution

anova(model_SoilCover) # Only Date has an effect 


# Post-hoc
Date_1<- subset(BB_data, BB_data$Date=="2018/06/07")
Date_2<- subset(BB_data, BB_data$Date=="2018/06/26")
Date_3<- subset(BB_data, BB_data$Date=="2018/07/17")
Date_4<- subset(BB_data, BB_data$Date=="2018/08/16")
Date_5<- subset(BB_data, BB_data$Date=="2019/09/06")

par(mfrow=c(1,5))
hist(Date_1$Soil_coverage_sqrt)
hist(Date_2$Soil_coverage_sqrt)
hist(Date_3$Soil_coverage_sqrt)
hist(Date_4$Soil_coverage_sqrt)
hist(Date_5$Soil_coverage_sqrt)

# Date 1 
model_SoilCover_D1 = lmer(Date_1$Soil_coverage_sqrt ~ Burning_intensity + (1|Block), data=Date_1,REML=TRUE)
plot(fitted(model_SoilCover_D1),residuals(model_SoilCover_D1))
anova(model_SoilCover_D1)

post_hoc_model_SoilCover_D1<-glht(model_SoilCover_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_SoilCover_D1)
tuk.cld1.SoilCover <- cld(post_hoc_model_SoilCover_D1)

# Date 2
model_SoilCover_D2 = lmer(Date_2$Soil_coverage_sqrt ~ Burning_intensity + (1|Block), data=Date_2,REML=TRUE)
plot(fitted(model_SoilCover_D2),residuals(model_SoilCover_D2))
anova(model_SoilCover_D2)

post_hoc_model_SoilCover_D2<-glht(model_SoilCover_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_SoilCover_D2)
tuk.cld2.SoilCover <- cld(post_hoc_model_SoilCover_D2)

# Date 3
model_SoilCover_D3 = lmer(Date_3$Soil_coverage_sqrt ~ Burning_intensity + (1|Block), data=Date_3,REML=TRUE)
plot(fitted(model_SoilCover_D3),residuals(model_SoilCover_D3))
anova(model_SoilCover_D3)

post_hoc_model_SoilCover_D3<-glht(model_SoilCover_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_SoilCover_D3)
tuk.cld3.SoilCover <- cld(post_hoc_model_SoilCover_D3)

# Date 4
model_SoilCover_D4 = lmer(Date_4$Soil_coverage_sqrt ~ Burning_intensity + (1|Block), data=Date_4,REML=TRUE)
plot(fitted(model_SoilCover_D4),residuals(model_SoilCover_D4))
anova(model_SoilCover_D4)

post_hoc_model_SoilCover_D4<-glht(model_SoilCover_D4, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_SoilCover_D4)
tuk.cld4.SoilCover <- cld(post_hoc_model_SoilCover_D4)

# Date 5
model_SoilCover_D5 = lmer(Date_5$Soil_coverage_sqrt ~ Burning_intensity + (1|Block), data=Date_5,REML=TRUE)
plot(fitted(model_SoilCover_D5),residuals(model_SoilCover_D5))
anova(model_SoilCover_D5)

post_hoc_model_SoilCover_D5<-glht(model_SoilCover_D5, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_SoilCover_D5)
tuk.cld5.SoilCover <- cld(post_hoc_model_SoilCover_D5)

Summarized_SoilCover <-BB_data %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(Soil_coverage, na.rm=TRUE), Sd = sd(Soil_coverage,na.rm=T),Max_Soil_coverage=max(Soil_coverage))

dat_text <- data.frame(label = c(tuk.cld1.SoilCover$mcletters$Letters, tuk.cld2.SoilCover$mcletters$Letters, tuk.cld3.SoilCover$mcletters$Letters, tuk.cld4.SoilCover$mcletters$Letters, tuk.cld5.SoilCover$mcletters$Letters),
                       Date   = c(rep("2018/06/07",4), rep("2018/06/26",4), rep("2018/07/17",4), rep("2018/08/16",4), rep("2019/09/06",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),5))

Summarized_SoilCov_tuk <- merge(Summarized_SoilCover,dat_text,by=c("Date","Burning_intensity"))


#### Soil coverage boxplot #### 

Date.labs <- c("23 days after burn","1 month after burn", "2 months after burn","3 months after burn","16 months after burn" )
names(Date.labs) <- c("2018/06/07","2018/06/26","2018/07/17","2018/08/16","2019/09/06")


Burning_Soil_coverage <- ggplot(BB_data, aes(x=Date,y=Soil_coverage))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Soil_coverage,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Soil_coverage,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_BBCov_tuk,aes(x=Date,y=2+Max_Soil_coverage,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  ylab("Soil coverage (%)")+
  guides(fill = "none")
Burning_Soil_coverage


##### Blueberry coverage #####

BB_data$BB_coverage_TR <- asin(sqrt((BB_data$BB_coverage/100)))

BB_hist_TR <- hist(BB_data$BB_coverage_TR)

model_BBCover = lmer(BB_coverage_TR ~ Burning_intensity * Date + (1|Block/Plot), data=BB_data,REML=TRUE)
plot(fitted(model_BBCover),residuals(model_BBCover)) # Needs to have no particular pattern

anova(model_BBCover) # Only Date has an effect 


# Post-hoc
Date_1<- subset(BB_data, BB_data$Date=="2018/06/07")
Date_2<- subset(BB_data, BB_data$Date=="2018/06/26")
Date_3<- subset(BB_data, BB_data$Date=="2018/07/17")
Date_4<- subset(BB_data, BB_data$Date=="2018/08/16")
Date_5<- subset(BB_data, BB_data$Date=="2019/09/06")

par(mfrow=c(1,5))
hist(Date_1$BB_coverage_TR)
hist(Date_2$BB_coverage_TR)
hist(Date_3$BB_coverage_TR)
hist(Date_4$BB_coverage_TR)
hist(Date_5$BB_coverage_TR)

# Date 1 
model_BBCover_D1 = lmer(Date_1$BB_coverage_TR ~ Burning_intensity + (1|Block), data=Date_1,REML=TRUE)
plot(fitted(model_BBCover_D1),residuals(model_BBCover_D1))
anova(model_BBCover_D1)

post_hoc_model_BBCover_D1<-glht(model_BBCover_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBCover_D1)
tuk.cld1.BBCover <- cld(post_hoc_model_BBCover_D1)

# Date 2
model_BBCover_D2 = lmer(Date_2$BB_coverage_TR ~ Burning_intensity + (1|Block), data=Date_2,REML=TRUE)
plot(fitted(model_BBCover_D2),residuals(model_BBCover_D2))
anova(model_BBCover_D2)

post_hoc_model_BBCover_D2<-glht(model_BBCover_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBCover_D2)
tuk.cld2.BBCover <- cld(post_hoc_model_BBCover_D2)

# Date 3
model_BBCover_D3 = lmer(Date_3$BB_coverage_TR ~ Burning_intensity + (1|Block), data=Date_3,REML=TRUE)
plot(fitted(model_BBCover_D3),residuals(model_BBCover_D3))
anova(model_BBCover_D3)

post_hoc_model_BBCover_D3<-glht(model_BBCover_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBCover_D3)
tuk.cld3.BBCover <- cld(post_hoc_model_BBCover_D3)

# Date 4
model_BBCover_D4 = lmer(Date_4$BB_coverage_TR ~ Burning_intensity + (1|Block), data=Date_4,REML=TRUE)
plot(fitted(model_BBCover_D4),residuals(model_BBCover_D4))
anova(model_BBCover_D4)

post_hoc_model_BBCover_D4<-glht(model_BBCover_D4, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBCover_D4)
tuk.cld4.BBCover <- cld(post_hoc_model_BBCover_D4)

# Date 5
model_BBCover_D5 = lmer(Date_5$BB_coverage_TR ~ Burning_intensity + (1|Block), data=Date_5,REML=TRUE)
plot(fitted(model_BBCover_D5),residuals(model_BBCover_D5))
anova(model_BBCover_D5)

post_hoc_model_BBCover_D5<-glht(model_BBCover_D5, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBCover_D5)
tuk.cld5.BBCover <- cld(post_hoc_model_BBCover_D5)

Summarized_BBCover <-BB_data %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(BB_coverage, na.rm=TRUE), Sd = sd(BB_coverage,na.rm=T),Max_BB_coverage=max(BB_coverage))

dat_text <- data.frame(label = c(tuk.cld1.BBCover$mcletters$Letters, tuk.cld2.BBCover$mcletters$Letters, tuk.cld3.BBCover$mcletters$Letters, tuk.cld4.BBCover$mcletters$Letters, tuk.cld5.BBCover$mcletters$Letters),
                       Date   = c(rep("2018/06/07",4), rep("2018/06/26",4), rep("2018/07/17",4), rep("2018/08/16",4), rep("2019/09/06",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),5))

Summarized_BBCov_tuk <- merge(Summarized_BBCover,dat_text,by=c("Date","Burning_intensity"))


#### Blueberry coverage boxplot #### 
Burning_BB_coverage <- ggplot(BB_data, aes(x=Date,y=BB_coverage))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=BB_coverage,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=BB_coverage,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_BBCov_tuk,aes(x=Date,y=2+Max_BB_coverage,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  ylab("Blueberry coverage (%)")+
  guides(fill = "none")
Burning_BB_coverage




##### Blueberry biomass #####

BB_data$BB_biomass_sqrt <- sqrt(BB_data$BB_biomass)

BBbiomass_hist_sqrt <- hist(BB_data$BB_biomass_sqrt)

model_BBbiomass= lmer(BB_biomass_sqrt ~ Burning_intensity * Date + (1|Block/Plot), data=BB_data,REML=TRUE)
plot(fitted(model_BBbiomass),residuals(model_BBbiomass)) # Needs to have no particular pattern

anova(model_BBbiomass) 

Date_1<- subset(BB_data, BB_data$Date=="2018/06/07")
Date_2<- subset(BB_data, BB_data$Date=="2018/06/26")
Date_3<- subset(BB_data, BB_data$Date=="2018/07/17")
Date_4<- subset(BB_data, BB_data$Date=="2018/08/16")
Date_5<- subset(BB_data, BB_data$Date=="2019/09/06")

par(mfrow=c(1,5))
hist(Date_1$BB_biomass_sqrt)
hist(Date_2$BB_biomass_sqrt)
hist(Date_3$BB_biomass_sqrt)
hist(Date_4$BB_biomass_sqrt)
hist(Date_5$BB_biomass_sqrt)
par(mfrow=c(1,1))

# Date 1 
model_BBbiomass_D1 = lmer(Date_1$BB_biomass_sqrt ~ Burning_intensity + (1|Block), data=Date_1,REML=TRUE)
plot(fitted(model_BBbiomass_D1),residuals(model_BBbiomass_D1))
anova(model_BBbiomass_D1)

post_hoc_model_BBbiomass_D1<-glht(model_BBbiomass_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBbiomass_D1)
tuk.cld1.BBbiomass <- cld(post_hoc_model_BBbiomass_D1)

# Date 2
model_BBbiomass_D2 = lmer(Date_2$BB_biomass_sqrt ~ Burning_intensity + (1|Block), data=Date_2,REML=TRUE)
plot(fitted(model_BBbiomass_D2),residuals(model_BBbiomass_D2))
anova(model_BBbiomass_D2)

post_hoc_model_BBbiomass_D2<-glht(model_BBbiomass_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBbiomass_D2)
tuk.cld2.BBbiomass <- cld(post_hoc_model_BBbiomass_D2)

# Date 3
model_BBbiomass_D3 = lmer(Date_3$BB_biomass_sqrt ~ Burning_intensity + (1|Block), data=Date_3,REML=TRUE)
plot(fitted(model_BBbiomass_D3),residuals(model_BBbiomass_D3))
anova(model_BBbiomass_D3)

post_hoc_model_BBbiomass_D3<-glht(model_BBbiomass_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBbiomass_D3)
tuk.cld3.BBbiomass <- cld(post_hoc_model_BBbiomass_D3)

# Date 4
model_BBbiomass_D4 = lmer(Date_4$BB_biomass_sqrt ~ Burning_intensity + (1|Block), data=Date_4,REML=TRUE)
plot(fitted(model_BBbiomass_D4),residuals(model_BBbiomass_D4))
anova(model_BBbiomass_D4)

post_hoc_model_BBbiomass_D4<-glht(model_BBbiomass_D4, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBbiomass_D4)
tuk.cld4.BBbiomass <- cld(post_hoc_model_BBbiomass_D4)

# Date 5
model_BBbiomass_D5 = lmer(Date_5$BB_biomass_sqrt ~ Burning_intensity + (1|Block), data=Date_5,REML=TRUE)
plot(fitted(model_BBbiomass_D5),residuals(model_BBbiomass_D5))
anova(model_BBbiomass_D5)

post_hoc_model_BBbiomass_D5<-glht(model_BBbiomass_D5, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBbiomass_D5)
tuk.cld5.BBbiomass <- cld(post_hoc_model_BBbiomass_D5)




Summarized_BBbiomass <-BB_data %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(BB_biomass, na.rm=TRUE), Sd = sd(BB_biomass,na.rm=T),Max_BB_biomass=max(BB_biomass))

dat_text <- data.frame(label = c(tuk.cld1.BBbiomass$mcletters$Letters, tuk.cld2.BBbiomass$mcletters$Letters, tuk.cld3.BBbiomass$mcletters$Letters, tuk.cld4.BBbiomass$mcletters$Letters, tuk.cld5.BBbiomass$mcletters$Letters),
                       Date   = c(rep("2018/06/07",4), rep("2018/06/26",4), rep("2018/07/17",4), rep("2018/08/16",4), rep("2019/09/06",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),5))

Summarized_BBbiomass_tuk <- merge(Summarized_BBbiomass,dat_text,by=c("Date","Burning_intensity"))




#### Blueberry biomass boxplot  #### 


Date.labs <- c("23 days after burn","1 month after burn", "2 months after burn","3 months after burn","16 months after burn" )
names(Date.labs) <- c("2018/06/07","2018/06/26","2018/07/17","2018/08/16","2019/09/06")

Burning_BB_biomass <- ggplot(BB_data, aes(x=Date,y=BB_biomass))+
                      geom_boxplot(aes(color=Burning_intensity))+
                      stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=BB_biomass,fill=Burning_intensity),
                                   position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
                      stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=BB_biomass,color=Burning_intensity),
                                   position=position_dodge(width=0.75))+
                      geom_text(data=Summarized_BBbiomass_tuk,aes(x=Date,y=10+Max_BB_biomass,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
                      facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
                      theme_light()+
                      theme(text =  element_text(size=15))+
                      theme(legend.position = "bottom",
                            legend.background = element_rect(size=0.4, 
                                                             linetype="solid", 
                                                             colour ="black"),
                            legend.text = element_text(size=15),
                            legend.title = element_text(size=16,face="bold"))+
                      theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
                      theme(axis.title.x = element_blank())+
                      scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
                      labs(y=expression(paste("Blueberry biomass ","(",g/m^2,")")))+
                      guides(fill = "none")
Burning_BB_biomass

####______________####
#### FIGURE 1 ####
ggarrange(Burning_BB_coverage, Burning_BB_biomass,
          ncol=1, nrow=2,
          labels="AUTO",font.label=list(size=12,face="bold"),
          common.legend = TRUE, legend="bottom")

