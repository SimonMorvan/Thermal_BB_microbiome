# 19_10_2021
# Simon Morvan



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 
#              Analysis of the evolution of several weeds (pin frame method)
#
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

### PACKAGES  ####
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(lme4)
library(lmerTest)
library(multcomp)

#### DATA ####


setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/Data/Weeds/")
Weeds_data<- read.csv2("Weed_evol_2018_2019.csv",dec=".")
Weeds_data$Burning_intensity<- factor(Weeds_data$Burning_intensity , levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))
Weeds_data$Date <- as.factor(Weeds_data$Date)
Weeds_data$Date<- factor(Weeds_data$Date, levels = c("2018/06/07","2018/06/26","2018/07/17","2018/08/16","2019/09/06"))


####______________####
###### RAW DATA DISTRIBUTION #####
####______________####

#### M. canadense coverage ####
Weeds_data$MC_coverage<- round(Weeds_data$MC_coverage,digits = 0)
Weeds_data$MC_coverage<- as.integer(Weeds_data$MC_coverage)
MC_hist <- hist(Weeds_data$MC_coverage)


#### G. procumbens coverage ####
# Take only the 3 last dates for GP
Weeds_dataGP <- subset(Weeds_data,Weeds_data$Date%in%c("2018/07/17","2018/08/16","2019/09/06"))
Weeds_dataGP$GP_coverage<- round(Weeds_dataGP$GP_coverage,digits = 0)
Weeds_dataGP$GP_coverage<- as.integer(Weeds_dataGP$GP_coverage)
GP_hist <- hist(Weeds_dataGP$GP_coverage)

mean_GP <-Weeds_data%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(GP_coverage, na.rm=TRUE))

#### C. canadensis ##### 
Weeds_dataCC<- subset(Weeds_data,Weeds_data$Date!="2018/06/07")
Weeds_dataCC$CC_coverage<- round(Weeds_dataCC$CC_coverage,digits = 0)
Weeds_dataCC$CC_coverage<- as.integer(Weeds_dataCC$CC_coverage)
CC_hist <- hist(Weeds_dataCC$CC_coverage)

mean_CC <-Weeds_data%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(CC_coverage, na.rm=TRUE))

### Total weed coverage ####
Weeds_data$TotWeed_coverage <- round(Weeds_data$TotWeed_coverage ,digits = 0)
Weeds_data$TotWeed_coverage<- as.integer(Weeds_data$TotWeed_coverage)
TotWeed_hist <- hist(Weeds_data$TotWeed_coverage)

par(mfrow=c(3,2))

plot(Soil_hist)
plot(Blueberry_hist)
plot(MC_hist)
plot(GP_hist)
plot(CC_hist)
plot(TotWeed_hist)




####______________####
###### MODELS  #####

####______________####

#### M. canadense #### 

# Poisson distribution --> GLMM 

GLMM_model_MC = glmer(MC_coverage ~ Burning_intensity * Date + (1|Block/Plot), 
                      data=Weeds_data, family="poisson")
summary(GLMM_model_MC)

plot(fitted(GLMM_model_MC),residuals(GLMM_model_MC)) # Needs to have no particular pattern


car::Anova(GLMM_model_MC,type=3) 


#### G. procumbens #### 

GLMM_model_GP = glmer(GP_coverage ~ Burning_intensity * Date + (1|Block/Plot), 
                      data=Weeds_dataGP, family="poisson")
summary(GLMM_model_GP)

plot(fitted(GLMM_model_GP),residuals(GLMM_model_GP)) # Needs to have no particular pattern

car::Anova(GLMM_model_GP,type=3) # No effect


#### C. canadensis  #### 

GLMM_model_CC = glmer(CC_coverage ~ Burning_intensity * Date + (1|Block/Plot), 
                      data=Weeds_dataCC, family="poisson")
summary(GLMM_model_CC)

plot(fitted(GLMM_model_CC),residuals(GLMM_model_CC)) # Needs to have no particular pattern

car::Anova(GLMM_model_CC,type=3) # Only data is significant 



### Total weed coverage ####
GLMM_model_TW = glmer(TotWeed_coverage ~ Burning_intensity * Date + (1|Block/Plot), 
                      data=Weeds_data, family="poisson")
summary(GLMM_model_TW)

plot(fitted(GLMM_model_TW),residuals(GLMM_model_TW)) # Needs to have no particular pattern

car::Anova(GLMM_model_TW,type=3) 

# Significant interaction between Burning intensity and date
# Significant date effect

# Normality of residuals 
hist(residuals(GLMM_model_TW)) # Histogram has to be bell shaped
qqnorm(residuals(GLMM_model_TW)) # Straight line = normal distribution



#### Total weed coverage per date analysis ##### 
Date_1<- subset(Weeds_data, Weeds_data$Date=="2018/06/07")
Date_2<- subset(Weeds_data, Weeds_data$Date=="2018/06/26")
Date_3<- subset(Weeds_data, Weeds_data$Date=="2018/07/17")
Date_4<- subset(Weeds_data, Weeds_data$Date=="2018/08/16")
Date_5<- subset(Weeds_data, Weeds_data$Date=="2019/09/06")


par(mfrow=c(1,5))
hist(Date_1$TotWeed_coverage)
hist(Date_2$TotWeed_coverage)
hist(Date_3$TotWeed_coverage)
hist(Date_4$TotWeed_coverage)
hist(Date_5$TotWeed_coverage)



# ANOVA per Date#
GLMM_TW_D1 <- glmer(TotWeed_coverage ~ Burning_intensity + (1|Block),data=Date_1, family="poisson")
plot(fitted(GLMM_TW_D1),residuals(GLMM_TW_D1)) # Needs to have no particular pattern

Anova_TW_D1 <- car::Anova(GLMM_TW_D1,type=3) # signif

post_hoc_TWD1<-glht(GLMM_TW_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_TWD1)
tuk.cld_TWD1 <- cld(post_hoc_TWD1)
# Significant difference between 3-1 and 3-2

GLMM_TW_D2 <- glmer(TotWeed_coverage ~ Burning_intensity + (1|Block),data=Date_2, family="poisson")
plot(fitted(GLMM_TW_D2),residuals(GLMM_TW_D2)) # Needs to have no particular pattern

car::Anova(GLMM_TW_D2,type=3) # Significant Burning Intensity !


post_hoc_TWD2<-glht(GLMM_TW_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_TWD2)
tuk.cld_TWD2 <- cld(post_hoc_TWD2)

## Difference is significant between 3-0 and 2_0 



GLMM_TW_D3 <- glmer(TotWeed_coverage ~ Burning_intensity + (1|Block),data=Date_3, family="poisson")
plot(fitted(GLMM_TW_D3),residuals(GLMM_TW_D3)) # Needs to have no particular pattern

car::Anova(GLMM_TW_D3,type=3) # Significant Burning Intensity !


post_hoc_TWD3<-glht(GLMM_TW_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_TWD3)
tuk.cld_TWD3 <- cld(post_hoc_TWD3)

## Difference is significant between 1-0, 2-0, 3-0, 2-1 and 3_1 



GLMM_TW_D4 <- glmer(TotWeed_coverage ~ Burning_intensity + (1|Block),data=Date_4, family="poisson")
plot(fitted(GLMM_TW_D4),residuals(GLMM_TW_D4)) # Needs to have no particular pattern

car::Anova(GLMM_TW_D4,type=3) # Significant Burning Intensity !


post_hoc_TWD4<-glht(GLMM_TW_D4, mcp(Burning_intensity="Tukey"))
summary(post_hoc_TWD4)
tuk.cld_TWD4 <- cld(post_hoc_TWD4)

## Difference is significant between 1-0, 2-0, 3-0, 2-1 and 3_1 


GLMM_TW_D5 <- glmer(TotWeed_coverage ~ Burning_intensity + (1|Block),data=Date_5, family="poisson")
plot(fitted(GLMM_TW_D5),residuals(GLMM_TW_D5)) # Needs to have no particular pattern

car::Anova(GLMM_TW_D5,type=3) # Significant Burning Intensity !


post_hoc_TWD5<-glht(GLMM_TW_D5, mcp(Burning_intensity="Tukey"))
summary(post_hoc_TWD5)
tuk.cld_TWD5 <- cld(post_hoc_TWD5)
## Difference is significant between 1-0, 2-0, 3-0, 2-1




####______________####
###### BOXPLOTS  #####
####______________####

#### M. canadense bp ####
Burning_MC <- ggplot(Weeds_data, aes(x=Date,y=MC_coverage))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=BB_coverage,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=BB_coverage,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  facet_grid(.~Date, drop=TRUE, scales="free_x",as.table = TRUE)+
  theme_light()+
  theme(strip.text.y = element_blank(),
        strip.background = element_blank())+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  ylab(expression(paste(italic("Maianthemum canadense"), " coverage (%)")))+
  guides(fill = "none")



#### G. procumbens coverage bp ####
Burning_GP<- ggplot(Weeds_dataGP, aes(x=Date,y=GP_coverage))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=BB_coverage,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=BB_coverage,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  facet_grid(.~Date, drop=TRUE, scales="free_x",as.table = TRUE)+
  theme_light()+
  theme(strip.text.y = element_blank(),
        strip.background = element_blank())+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  ylab(expression(paste(italic("Gaultheria procumbens"), " coverage (%)")))+
  guides(fill = "none")


#### C. canadensis  bp##### 
Weeds_dataCC$Date
Date.labs <- c("1 month after burn", "2 months after burn","3 months after burn","13 months after burn" )
names(Date.labs) <- c("2018/06/26","2018/07/17", "2018/08/16","2019/09/06")

Burning_CC <- ggplot(Weeds_dataCC, aes(x=Date,y=CC_coverage))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=CC_coverage,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=CC_coverage,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  theme_light()+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(text =  element_text(size=14))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  ylab(expression(paste(italic("Cornus canadensis"), " coverage (%)")))+
  guides(fill = "none")
Burning_CC


Burning_CC_point <- ggplot(Weeds_dataCC, aes(x=Date,y=CC_coverage))+
  geom_jitter(aes(color=Burning_intensity),position = position_jitter(width=0.3,seed = 0.2)) +
  geom_text_repel(aes(label = Plot),size=4,max.overlaps = 20,position = position_jitter(width=0.3,seed = 0.2))+
  theme(text =  element_text(size=15))+
  facet_grid(.~Date, drop=TRUE, scales="free_x",as.table = TRUE)+
  theme_light()+
  theme(strip.text.y = element_blank(),
        strip.background = element_blank())+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  ylab(expression(paste(italic("Cornus canadensis"), " coverage (%)")))+
  guides(fill = "none")
Burning_CC_point



### Total weed coverage bp ####

Summarized_TW <-Weeds_data %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(TotWeed_coverage, na.rm=TRUE), Sd = sd(TotWeed_coverage,na.rm=T),Max_TotWeed_coverage=max(TotWeed_coverage))

dat_text <- data.frame(label = c(tuk.cld_TWD1$mcletters$Letters, tuk.cld_TWD2$mcletters$Letters, tuk.cld_TWD3$mcletters$Letters, tuk.cld_TWD4$mcletters$Letters, tuk.cld_TWD5$mcletters$Letters),
                       Date   = c(rep("2018/06/07",4), rep("2018/06/26",4), rep("2018/07/17",4), rep("2018/08/16",4), rep("2019/09/06",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),5))

Summarized_TW_tuk <- merge(Summarized_TW,dat_text,by=c("Date","Burning_intensity"))

Date.labs <- c("23 days after burn","1 month after burn", "2 months after burn","3 months after burn","16 months after burn" )
names(Date.labs) <- c("2018/06/07","2018/06/26","2018/07/17", "2018/08/16","2019/09/06")


Burning_TotWeed <- ggplot(Weeds_data, aes(x=Date,y=TotWeed_coverage))+
  geom_boxplot(aes(color=Burning_intensity))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=TotWeed_coverage,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=TotWeed_coverage,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_TW_tuk,aes(x=Date,y=2+Max_TotWeed_coverage,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  theme_light()+
  theme(text =  element_text(size=14))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  scale_colour_manual(name= "Burning intensity", values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  ylab("Total weeds coverage (%)")+
  guides(fill = "none")
Burning_TotWeed


####______________####
### FIGURE 2 #####
ggarrange(Burning_CC, Burning_TotWeed,
          ncol=1, nrow=2,
          labels="AUTO",font.label=list(size=12,face="bold"),
          common.legend = TRUE, legend="bottom")







