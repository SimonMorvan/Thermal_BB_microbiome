# 05_02_2021
# Simon Morvan

R.version #4.1.1 
rm(list = ls())
gc()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#
#             Analysis of the  blueberry                                                                     
#                                                                            
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

##### _____________ #####
#### PACKAGES ####
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lmerTest)
library(MASS)
library(multcomp)


##### _____________ #####
#### DATA ####
setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/Data/Bleuet/")
BB_2019 <- read.csv2("Blueberry_harvest.csv")
colnames(BB_2019)

BB_2019$Plot <- as.factor(as.character(BB_2019$Plot))
BB_2019$Burning_intensity<- factor(BB_2019$Burning_intensity , levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))
BB_2019$Block <- as.factor(as.character(BB_2019$Block))

##### _____________ #####
####  SHOOTS #### 
#### Shoot number ####
Burning_bb_Shoot_number <- ggplot(BB_2019, aes(y=Shoot_Total_Number_.1m2., x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Shoot_Total_Number_.1m2.,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Shoot_Total_Number_.1m2.,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Shoot number per square meter")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")
Burning_bb_Shoot_number

ShootnumberHist <- hist(BB_2019$Shoot_Total_Number_.1m2.)

mean_shootdensity <-BB_2019 %>%
  group_by(Burning_intensity) %>%
  summarize(Mean = mean(Shoot_Total_Number_.1m2., na.rm=TRUE))


Model_ShootNumber<- lmer(Shoot_Total_Number_.1m2. ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_ShootNumber),residuals(Model_ShootNumber))
anova(Model_ShootNumber)
# Non significant


##### Shoot height #####
Burning_bb_Shoot_height <- ggplot(BB_2019, aes(y=Mean_height_9_shoots_.cm., x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Mean_height_9_shoots_.cm.,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Mean_height_9_shoots_.cm.,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Mean shoot height (cm)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_bb_Shoot_height

ShootHeightHist <- hist(BB_2019$Mean_height_9_shoots_.cm.)

mean_shootheight <-BB_2019 %>%
  group_by(Burning_intensity) %>%
  summarize(Mean = mean(Mean_height_9_shoots_.cm., na.rm=TRUE))


Model_ShootHeight <- lmer(Mean_height_9_shoots_.cm. ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_ShootHeight),residuals(Model_ShootHeight)) # Needs to have no particular pattern
anova(Model_ShootHeight)


Shoots <- ggarrange(Burning_bb_Shoot_number,
                    Burning_bb_Shoot_height,
                    ncol=2, 
                    nrow=2, 
                    common.legend = TRUE, 
                    legend="bottom")

par(mfrow=c(2,1))
plot(ShootNumberHist)
plot(ShootHeightHist)

##### _____________ #####
#### BLUEBERRIES #####

### Fresh weight ####
FreshFruitWeightHist <- hist(BB_2019$Fresh_fruit_weight_.T.ha.)

Model_FreshFruitWeight<- lmer(Fresh_fruit_weight_.T.ha. ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_FreshFruitWeight),residuals(Model_FreshFruitWeight)) # Needs to have no particular pattern
anova(Model_FreshFruitWeight)

post_hoc_model_BBYield<-glht(Model_FreshFruitWeight, mcp(Burning_intensity="Tukey"))
summary(post_hoc_model_BBYield)
tuk.cld.BBYields <- cld(post_hoc_model_BBYield)

yield <-BB_2019 %>%
  group_by(Burning_intensity) %>%
  summarize(Mean = mean(Fresh_fruit_weight_.T.ha.), Sd = sd(Fresh_fruit_weight_.T.ha.),Max_BBYields=max(Fresh_fruit_weight_.T.ha.))

dat_text <- data.frame(label = tuk.cld.BBYields$mcletters$Letters,
                       Burning_intensity  = c("NoBurn","LowBurn","MedBurn","HighBurn"))

yield_tuk <- merge(yield,dat_text,by="Burning_intensity")



Burning_bb_FW <- ggplot(BB_2019, aes(y=Fresh_fruit_weight_.T.ha., x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Fresh_fruit_weight_.T.ha.,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Fresh_fruit_weight_.T.ha.,color=Burning_intensity),
               position=position_dodge(width=0.75))+  
  geom_text(data=yield_tuk,aes(y=0.1+Max_BBYields,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
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
        legend.title = element_text(size=16,face="bold"))+
  labs(y=expression(paste("Blueberry fresh weight (",T.h^-1,")")))+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")
Burning_bb_FW



#### Blueberry number in 125ml #####
Burning_bb_number <- ggplot(BB_2019, aes(y=Total_fruit_number_.125ml._Blue_Green, x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Total_fruit_number_.125ml._Blue_Green,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Total_fruit_number_.125ml._Blue_Green,color=Burning_intensity),
               position=position_dodge(width=0.75))+  
  geom_text(data=,aes(x=Date,y=10+Max_BB_biomass,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Blueberry number in 125mL")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_bb_number

TotFruitNumHist <-hist(BB_2019$Total_fruit_number_.125ml._Blue_Green)

Model_TotFruitNum<- lmer(Total_fruit_number_.125ml._Blue_Green ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_TotFruitNum),residuals(Model_TotFruitNum)) # Needs to have no particular pattern
anova(Model_TotFruitNum)


Blueberry <- ggarrange(Burning_bb_FW,
                       Burning_bb_number,
                       ncol=1, 
                       nrow=2, 
                       common.legend = TRUE, 
                       legend="bottom")

par(mfrow=c(2,1))
plot(FreshFruitWeightHist)
plot(TotFruitNumHist)


##### Frozen blueberries ######
Burning_bb_frozen_avweight <- ggplot(BB_2019, aes(y=Average_weight_frozen_fruit_Blue_Green_.g., x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Average_weight_frozen_fruit_Blue_Green_.g.,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Average_weight_frozen_fruit_Blue_Green_.g.,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Frozen blueberry average weight (g)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_bb_frozen_avweight


hist(BB_2019$Average_weight_frozen_fruit_Blue_Green_.g.)

Model_FrozenAvWeight<- lmer(Average_weight_frozen_fruit_Blue_Green_.g. ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_FrozenAvWeight),residuals(Model_FrozenAvWeight)) # Needs to have no particular pattern
anova(Model_FrozenAvWeight) #Non significant

mean_frozenweight <-BB_2019 %>%
                    group_by(Burning_intensity) %>%
                    summarize(Mean = mean(Average_weight_frozen_fruit_Blue_Green_.g., na.rm=TRUE))

##### Dried blueberries ######
Burning_bb_avDW <- ggplot(BB_2019, aes(y=Average_weight_dried_fruit_Blue_Green_.g., x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Average_weight_dried_fruit_Blue_Green_.g.,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Average_weight_dried_fruit_Blue_Green_.g.,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Dried blueberry average weight (g)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_bb_avDW

hist(BB_2019$Dried_fruit_weight_.g.125ml._Blue_Green)

Model_Dryweight<- lmer(Average_weight_dried_fruit_Blue_Green_.g. ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_Dryweight),residuals(Model_Dryweight))
anova(Model_Dryweight)
# Non significant






##### _____________ #####
##### BLUE BLUEBERRIES######

### Ripe Blueberry %  ####
BB_2019$RipeProportion <- BB_2019$Blue_fruits_number/BB_2019$Total_fruit_number_.125ml._Blue_Green

Burning_RipeProportion<- ggplot(BB_2019, aes(y=RipeProportion, x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=RipeProportion,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=RipeProportion,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Ripe blueberry proportion")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_RipeProportion

mean_ripeprop <-BB_2019 %>%
  group_by(Burning_intensity) %>%
  summarize(Mean = mean(RipeProportion, na.rm=TRUE))

hist(BB_2019$RipeProportion)

Model_RipeProp<- lmer(RipeProportion ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_RipeProp),residuals(Model_RipeProp))
anova(Model_RipeProp)




### Blue Blueberry number in 125ml ####
Burning_Blue_bb_number <- ggplot(BB_2019, aes(y=Blue_fruits_number, x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Blue_fruits_number,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Blue_fruits_number,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Blue blueberry number in 125mL")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_Blue_bb_number

hist(BB_2019$Blue_fruits_number)

Model_BlueNum<- lmer(Blue_fruits_number ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_BlueNum),residuals(Model_BlueNum))
anova(Model_BlueNum)


##### Blue Blueberry average weigh #####
Burning_Blue_bb_avweigh <- ggplot(BB_2019, aes(y=Average_blue_fruits_weight_.g., x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Average_blue_fruits_weight_.g.,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Average_blue_fruits_weight_.g.,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Blue blueberry average weigh (g)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")
Burning_Blue_bb_avweigh

mean_ripeweight <-BB_2019 %>%
           group_by(Burning_intensity) %>%
           summarize(Mean = mean(Average_blue_fruits_weight_.g., na.rm=TRUE))



BlueBBWeigthHist <- hist(BB_2019$Average_blue_fruits_weight_.g.)

Model_bBBweigh <- lmer(Average_blue_fruits_weight_.g. ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_bBBweigh),residuals(Model_bBBweigh)) # Needs to have no particular pattern
anova(Model_bBBweigh)

BlueBBWeigthHistLog <- hist(log(BB_2019$Average_blue_fruits_weight_.g.))
Model_bBBweighLog <- lmer(log(Average_blue_fruits_weight_.g.) ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_bBBweighLog),residuals(Model_bBBweighLog)) # Needs to have no particular pattern
anova(Model_bBBweighLog)

BlueBBWeigthHistSqrt<- hist(sqrt(BB_2019$Average_blue_fruits_weight_.g.))
Model_bBBweighSqrt <- lmer(sqrt(Average_blue_fruits_weight_.g.) ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_bBBweighSqrt),residuals(Model_bBBweighSqrt)) # Needs to have no particular pattern
anova(Model_bBBweighSqrt)



bc <- boxcox(BB_2019$Average_blue_fruits_weight_.g.~ BB_2019$Burning_intensity)
(lambda <- bc$x[which.max(bc$y)])

Model_bBBweighBoxCox <- lmer(((Average_blue_fruits_weight_.g.^lambda-1)/lambda) ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_bBBweighBoxCox),residuals(Model_bBBweighBoxCox)) # Needs to have no particular pattern
anova(Model_bBBweighBoxCox)


post_hoc_RipeWeight<-glht(Model_bBBweighBoxCox, mcp(Burning_intensity="Tukey"))
summary(post_hoc_RipeWeight)
tuk.cld <- cld(post_hoc_RipeWeight)


Blue_Blueberry <- ggarrange(Burning_Blue_bb_number,
                            Burning_Blue_bb_avweigh,
                            ncol=1, 
                            nrow=2, 
                            common.legend = TRUE, 
                            legend="bottom")

par(mfrow=c(2,1))
plot(BlueBBNumHist)
plot(BlueBBWeigthHist)


##### _____________ #####
##### GREEN BLUEBERRIES######

#### Green Blueberry number in 125ml #####
Burning_Green_bb_number <- ggplot(BB_2019, aes(y=Green_fruits_number, x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Green_fruits_number,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Green_fruits_number,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Green blueberry number in 125mL")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_Green_bb_number

GreenBBNumHist <- hist(BB_2019$Green_fruits_number)

Model_GreenNum<- lmer(Green_fruits_number ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_GreenNum),residuals(Model_GreenNum))
anova(Model_GreenNum)


##### Green Blueberry average weight #####
Burning_Green_bb_avweight <- ggplot(BB_2019, aes(y=Average_green_fruits_weight_.g., x=Burning_intensity)) + 
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(y=Average_green_fruits_weight_.g.,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(y=Average_green_fruits_weight_.g.,color=Burning_intensity),
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
        legend.title = element_text(size=16,face="bold"))+
  ylab("Green blueberry average weigh (g)")+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  guides(fill = "none")

Burning_Green_bb_avweight

GreenBBWeigthHist <- hist(BB_2019$Average_green_fruits_weight_.g.)

Model_GreenWeight<- lmer(Average_green_fruits_weight_.g. ~ Burning_intensity + (1|Block), data=BB_2019,REML=TRUE)
plot(fitted(Model_GreenWeight),residuals(Model_GreenWeight))
anova(Model_GreenWeight)

Green_Blueberry <- ggarrange(Burning_Green_bb_number,
                             Burning_Green_bb_avweight,
                             ncol=1, 
                             nrow=2, 
                             common.legend = TRUE, 
                             legend="bottom")



par(mfrow=c(2,1))
plot(GreenBBNumHist)
plot(GreenBBWeigthHist)
