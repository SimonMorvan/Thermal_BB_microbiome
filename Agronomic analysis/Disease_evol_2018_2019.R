# 01_02_2021
# Simon Morvan
# 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 
# Analysis of the evolution of several disease (pin frame method)
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

### PACKAGES 
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lmerTest)
library(multcomp)

## DATA #####
setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/Data/Disease/")
Disease<- read.csv2("Disease_evol_2018_2019.csv",dec=".")
colnames(Disease)


Disease$Plot <- as.factor(as.character(Disease$Plot))
Disease$Burning_intensity<- factor(Disease$Burning_intensity , levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))
Disease$Block <- as.factor(as.character(Disease$Block))
Disease$Date<- factor(Disease$Date, levels = c("2018/07/17","2018/08/16","2019/09/06"))


Date_1<- subset(Disease, Disease$Date=="2018/07/17")
Date_2<- subset(Disease, Disease$Date=="2018/08/16")
Date_3<- subset(Disease, Disease$Date=="2019/09/06")


##### _______ #####
###### Healthy #####
Healthy_data <- subset(Disease,Disease$Date%in%c("2018/07/17","2018/08/16"))

Burning_Healthy <- ggplot(Healthy_data, aes(x=Date,y=Healthy))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Healthy,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Healthy,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  ylab("Healthy plants (%)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")
  
Burning_Healthy


hist(Disease$Healthy)
Disease$HealthyTR <- asin(sqrt((Disease$Healthy/100)))
hist(Disease$HealthyTR)

model_Healthy = lmer(HealthyTR ~ Burning_intensity * Date + (1|Block/Plot), data=Disease,REML=TRUE)

plot(fitted(model_Healthy),residuals(model_Healthy)) # Needs to have no particular pattern
anova(model_Healthy) 


hist(residuals(model_Healthy))
# Histogram has to be bell shaped
qqnorm(residuals(model_Healthy))
# Straight line = normal distribution



##### _______ #####
###### Red tip on bottom leaf #####

Burning_RTBL <- ggplot(Disease, aes(x=Date,y=Red_tip_bottom_leaf))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Red_tip_bottom_leaf,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Red_tip_bottom_leaf,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  ylab("Red tip on bottom leaf (%)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_RTBL

RTBL_data <- subset(Disease,Disease$Date%in%c("2018/07/17","2018/08/16"))

hist(RTBL_data$Red_tip_bottom_leaf)
RTBL_data$Red_tip_bottom_leaf <- round(RTBL_data$Red_tip_bottom_leaf,digits = 0)
RTBL_data$Red_tip_bottom_leaf <- as.integer(RTBL_data$Red_tip_bottom_leaf)
hist(RTBL_data$Red_tip_bottom_leaf)


GLMM_model_RTBL= glmer(Red_tip_bottom_leaf ~ Burning_intensity * Date + (1|Block/Plot), 
                      data=RTBL_data, family="poisson")
summary(GLMM_model_RTBL)

plot(fitted(GLMM_model_RTBL),residuals(GLMM_model_RTBL)) # Needs to have no particular pattern


car::Anova(GLMM_model_RTBL,type=3) # Only data is significant
# Burning intensity pval= 0.055 / Interaction = 0.06

###### Date per date #####

RTBL_date1 <- subset(RTBL_data,RTBL_data$Date%in%c("2018/07/17"))
RTBL_date2 <- subset(RTBL_data,RTBL_data$Date%in%c("2018/08/16"))

par(mfrow=c(1,2))
hist(RTBL_date1$Red_tip_bottom_leaf)
hist(RTBL_date2$Red_tip_bottom_leaf)

# ANOVA per Date#
GLMM_RTBL_D1 <- glmer(Red_tip_bottom_leaf ~ Burning_intensity + (1|Block),data=RTBL_date1, family="poisson")

plot(fitted(GLMM_RTBL_D1),residuals(GLMM_RTBL_D1)) # Needs to have no particular pattern

Anova_RTBL_D1 <- car::Anova(GLMM_RTBL_D1,type=3) # signif

post_hoc_RTBLD1<-glht(GLMM_RTBL_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_RTBLD1)
tuk.cld1 <- cld(post_hoc_RTBLD1)

# Significant difference between 1-0 and 3-0; 2-1 ; 3-2 


GLMM_RTBL_D2 <- glmer(Red_tip_bottom_leaf ~ Burning_intensity + (1|Block),data=RTBL_date2, family="poisson")

plot(fitted(GLMM_RTBL_D2),residuals(GLMM_RTBL_D2)) # Needs to have no particular pattern

Anova_RTBL_D2 <- car::Anova(GLMM_RTBL_D2,type=3) # NON-signif

post_hoc_RTBLD2<-glht(GLMM_RTBL_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_RTBLD2)



##### _______ #####
###### Black leaf bud NOT ENOUGH DATA #####

Burning_BLB <- ggplot(Disease, aes(x=Date,y=Black_leaf_bud,color=Burning_intensity))+
  geom_boxplot()+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  ylab("Black leaf bud (%)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))
Burning_BLB


BLB_data <- subset(Disease,Disease$Date%in%c("2018/07/17","2018/08/16"))

hist(BLB_data$Black_leaf_bud)
BLB_data$Black_leaf_bud <- round(BLB_data$Black_leaf_bud,digits = 0)
BLB_data$Black_leaf_bud <- as.integer(BLB_data$Black_leaf_bud)
hist(BLB_data$Black_leaf_bud)


GLMM_model_BLB= glmer(Black_leaf_bud ~ Burning_intensity * Date + (1|Block/Plot), 
                       data=BLB_data, family="poisson")
summary(GLMM_model_BLB)

plot(fitted(GLMM_model_BLB),residuals(GLMM_model_BLB)) # Needs to have no particular pattern


car::Anova(GLMM_model_RTBL,type=3) 


##### _______ #####
###### Dry leaf apex  #####

Burning_DLA <- ggplot(Disease, aes(x=Date,y=Dry_leaf_apex))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Dry_leaf_apex,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Dry_leaf_apex,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  ylab("Dry leaf apex (%)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_DLA

DLA_data <- subset(Disease,Disease$Date%in%c("2018/07/17","2018/08/16"))

hist(DLA_data$Dry_leaf_apex)

DLA_data$Dry_leaf_apex <- round(DLA_data$Dry_leaf_apex,digits = 0)
DLA_data$Dry_leaf_apex <- as.integer(DLA_data$Dry_leaf_apex)
hist(DLA_data$Dry_leaf_apex)

GLMM_model_DLA = glmer(Dry_leaf_apex ~ Burning_intensity * Date + (1|Block/Plot), 
                      data=DLA_data, family="poisson")

plot(fitted(GLMM_model_DLA),residuals(GLMM_model_DLA)) # Needs to have no particular pattern

car::Anova(GLMM_model_DLA,type=3)  



##### _______ #####
###### Septoria #####

Burning_Septoria <- ggplot(Disease, aes(x=Date,y=Septoria))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Septoria,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Septoria,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  ylab("Septoria (%)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_Septoria


hist(Disease$Septoria)
Disease$SeptoriaTR <- asin(sqrt((Disease$Septoria/100)))
hist(Disease$SeptoriaTR)


model_Septoria = lmer(SeptoriaTR ~ Burning_intensity * Date + (1|Block/Plot), data=Disease,REML=TRUE)

plot(fitted(model_Septoria),residuals(model_Septoria)) # Needs to have no particular pattern
# I. Linearity assumption 
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# III. Homoskedasticity
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Septoria))
# Histogram has to be bell shaped
qqnorm(residuals(model_Septoria))
# Straight line = normal distribution

anova(model_Septoria) # Only Date has an effect 


#### PER DATE ####
library(multcomp)

par(mfrow=c(1,3))
hist(Date_1$Septoria)
hist(Date_2$Septoria)
hist(Date_3$Septoria)

model_Septoria_D1 = lmer(Septoria ~ Burning_intensity + (1|Block), data=Date_1,REML=TRUE)
plot(fitted(model_Septoria_D1),residuals(model_Septoria_D1)) # Needs to have no particular pattern
anova(model_Septoria_D1)

post_hoc_model_Septoria_D1<-glht(model_Septoria_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_Septoria_D1)
tuk.cld1 <- cld(post_hoc_model_Septoria_D1)

model_Septoria_D2 = lmer(Septoria ~ Burning_intensity + (1|Block), data=Date_2,REML=TRUE)
plot(fitted(model_Septoria_D2),residuals(model_Septoria_D2)) # Needs to have no particular pattern
anova(model_Septoria_D2)

post_hoc_model_Septoria_D2<-glht(model_Septoria_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_Septoria_D2)
tuk.cld2 <- cld(post_hoc_model_Septoria_D2)


model_Septoria_D3 = lmer(Septoria ~ Burning_intensity + (1|Block), data=Date_3,REML=TRUE)
plot(fitted(model_Septoria_D3),residuals(model_Septoria_D3)) # Needs to have no particular pattern
anova(model_Septoria_D3)

post_hoc_model_Septoria_D3<-glht(model_Septoria_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_Septoria_D3)
tuk.cld3 <- cld(post_hoc_model_Septoria_D3)

Summarized_Septoria <-Disease %>%
                      group_by(Date,Burning_intensity) %>%
                      summarize(Mean = mean(Septoria, na.rm=TRUE), Sd = sd(Septoria,na.rm=T),Max_Septoria=max(Septoria))

dat_text <- data.frame(label = c(tuk.cld1$mcletters$Letters, tuk.cld2$mcletters$Letters,  tuk.cld3$mcletters$Letters),
                       Date   = c(rep("2018/07/17",4), rep("2018/08/16",4), rep("2019/09/06",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))

Summarized_Septoria_tuk <- merge(Summarized_Septoria,dat_text,by=c("Date","Burning_intensity"))



Date.labs <- c("2 months after burn", "3 months after burn","16 months after burn" )
names(Date.labs) <- c("2018/07/17","2018/08/16", "2019/09/06")

#### FIGURE 3 ######
Burning_Septoria <- ggplot(Disease, aes(x=Date,y=Septoria))+
  geom_boxplot(aes(color=Burning_intensity))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Septoria,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Septoria,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_Septoria_tuk,aes(x=Date,y=2+Max_Septoria,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
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
  ylab("Septoria leaf spot disease (%)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_Septoria



##### _______ #####
###### Red plant LACK OF DATA #####

Burning_RP <- ggplot(Disease, aes(x=Date,y=Red_plant))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Red_plant,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Red_plant,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  ylab("Red plant (%)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_RP

RedPlant_data <- subset(Disease,Disease$Date%in%c("2018/08/16"))

hist(RedPlant_data$Red_plant)

RedPlant_data$Red_plant <- round(RedPlant_data$Red_plant,digits = 0)
RedPlant_data$Red_plant <- as.integer(RedPlant_data$Red_plant)
hist(RedPlant_data$Red_plant)

GLMM_model_RP = glmer(Red_plant ~ Burning_intensity + (1|Block), 
                       data=RedPlant_data, family="poisson")

plot(fitted(GLMM_model_RP),residuals(GLMM_model_RP)) # Needs to have no particular pattern

car::Anova(GLMM_model_RP,type=3) # Only data is significant 


##### _______ #####
##### DISEASE PLOTS #####

Disease_Plots <- ggarrange(Burning_Healthy_mean,
                        Burning_Septoria_mean,
                        Burning_RTBL_mean,
                        Burning_RP_mean,
                        Burning_DLA_mean,
                        Burning_BLB_mean,
                        ncol=3, 
                        nrow=2, 
                        common.legend = TRUE, 
                        legend="bottom")


