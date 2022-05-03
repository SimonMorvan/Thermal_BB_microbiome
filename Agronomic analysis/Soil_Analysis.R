# 29_04_2022
# Simon Morvan
# 
##### Loading and formating ######
setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/Scripts_Github/Agricultural data/Soil/")
list.files()
Soil_chem<- read.csv2("Chemistry.csv", dec=".")
colnames(Soil_chem)

library(ggplot2)
library(ggpubr)
library(lmerTest)
library(multcomp)
library(dplyr)

Soil_chem$Plot <- as.factor(as.character(Soil_chem$Plot))
Soil_chem$Soil_depth <- as.factor(as.character(Soil_chem$Soil_depth))
Soil_chem$Burning_intensity<- factor(Soil_chem$Burning_intensity , levels = c("NoBurn", "LowBurn", "MedBurn","HighBurn"))
Soil_chem$Block <- as.factor(as.character(Soil_chem$Block))
Soil_chem$Date<- as.factor(as.character(Soil_chem$Date)) 

Soil_chem_Org <- subset(Soil_chem,Soil_chem$Soil_depth=="O")
Soil_chem_Min <- subset(Soil_chem,Soil_chem$Soil_depth=="M")

##### Separation per date ####
Soil_chem_Org_D1 <- subset(Soil_chem_Org,Soil_chem_Org$Date=="2018/06/15")
Soil_chem_Org_D2 <- subset(Soil_chem_Org,Soil_chem_Org$Date=="2018/09/15")
Soil_chem_Org_D3 <- subset(Soil_chem_Org,Soil_chem_Org$Date=="2019/09/24")

Soil_chem_Min_D1 <- subset(Soil_chem_Min,Soil_chem_Min$Date=="2018/06/15")
Soil_chem_Min_D2 <- subset(Soil_chem_Min,Soil_chem_Min$Date=="2018/09/15")
Soil_chem_Min_D3 <- subset(Soil_chem_Min,Soil_chem_Min$Date=="2019/09/24")


#####_______#####

###### pH Organic MM ####

hist(Soil_chem_Org$pH)

model_phOrg = lmer(pH ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Org,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_phOrg),residuals(model_phOrg))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_phOrg),residuals(model_phOrg))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_phOrg))
# Histogram has to be bell shaped
qqnorm(residuals(model_phOrg))
# Straight line = normal distribution


anova(model_phOrg)

# Test per date
# Date 1
model_ph_Org_D1 = lmer(pH ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_ph_Org_D1),residuals(model_ph_Org_D1)) # Needs to have no particular pattern
anova(model_ph_Org_D1)

post_hoc_ph_org_D1<-glht(model_ph_Org_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_ph_org_D1)
tuk.cld1.ph_org <- cld(post_hoc_ph_org_D1)

# Date 2
model_ph_Org_D2 = lmer(pH ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_ph_Org_D2),residuals(model_ph_Org_D2)) # Needs to have no particular pattern
anova(model_ph_Org_D2)

post_hoc_ph_org_D2<-glht(model_ph_Org_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_ph_org_D2)
tuk.cld2.ph_org <- cld(post_hoc_ph_org_D2)

# Date 3
model_ph_Org_D3 = lmer(pH ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_ph_Org_D3),residuals(model_ph_Org_D3)) # Needs to have no particular pattern
anova(model_ph_Org_D3)

post_hoc_ph_org_D3<-glht(model_ph_Org_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_ph_org_D3)
tuk.cld3.ph_org <- cld(post_hoc_ph_org_D3)


mean_pHorg<-Soil_chem_Org%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(pH, na.rm=TRUE))

###### pH Organic Boxplot ####

Date.labs <- c("1 month after burn","4 months after burn","16 months after burn")
names(Date.labs) <- c("2018/06/15", "2018/09/15","2019/09/24")

Burning_pH_evolution_Org <- ggplot(Soil_chem_Org, aes(x=Date,y=pH))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(3.8,5.4)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=pH,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=pH,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  ylab("pH")+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  guides(fill = "none")

Burning_pH_evolution_Org

###### pH Mineral MM ####

hist(Soil_chem_Min$pH)

model_phMin = lmer(pH ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Min,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_phMin),residuals(model_phMin))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_phMin),residuals(model_phMin))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_phMin))
# Histogram has to be bell shaped
qqnorm(residuals(model_phMin))
# Straight line = normal distribution

anova(model_phMin)

# Test per date
# Date 1
model_ph_Min_D1 = lmer(pH ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D1,REML=TRUE)
plot(fitted(model_ph_Min_D1),residuals(model_ph_Min_D1)) # Needs to have no particular pattern
anova(model_ph_Min_D1)

post_hoc_ph_Min_D1<-glht(model_ph_Min_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_ph_Min_D1)
tuk.cld1.ph_Min <- cld(post_hoc_ph_Min_D1)


# Date 2
model_ph_Min_D2 = lmer(pH ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D2,REML=TRUE)
plot(fitted(model_ph_Min_D2),residuals(model_ph_Min_D2)) # Needs to have no particular pattern
anova(model_ph_Min_D2)

post_hoc_ph_Min_D2<-glht(model_ph_Min_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_ph_Min_D2)
tuk.cld2.ph_Min <- cld(post_hoc_ph_Min_D2)


# Date 3
model_ph_Min_D3 = lmer(pH ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D3,REML=TRUE)
plot(fitted(model_ph_Min_D3),residuals(model_ph_Min_D3)) # Needs to have no particular pattern
anova(model_ph_Min_D3)

post_hoc_ph_Min_D3<-glht(model_ph_Min_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_ph_Min_D3)
tuk.cld3.ph_Min <- cld(post_hoc_ph_Min_D3)


mean_pHmin<-Soil_chem_Min%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(pH, na.rm=TRUE))

###### pH Mineral Boxplot ####

Burning_pH_evolution_Min <- ggplot(Soil_chem_Min, aes(x=Date,y=pH))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(3.8,5.4)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=pH,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=pH,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  ylab("pH")+
  guides(fill = "none")

Burning_pH_evolution_Min

#### Common plot ####
Burning_pH_evolution <- ggarrange(Burning_pH_evolution_Org,
                                  Burning_pH_evolution_Min,
                                  labels= c("A","B"),
                                  ncol=1, 
                                  nrow=2, 
                                  common.legend = TRUE, 
                                  legend="bottom")


#####_______#####

###### C Organic MM ####

hist(Soil_chem_Org$Carbon_PC)

model_Corg= lmer(Carbon_PC ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Org,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Corg),residuals(model_Corg))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Corg),residuals(model_Corg))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Corg))
# Histogram has to be bell shaped
qqnorm(residuals(model_Corg))
# Straight line = normal distribution


anova(model_Corg)


# Date 1
model_C_Org_D1 = lmer(Carbon_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_C_Org_D1),residuals(model_C_Org_D1)) # Needs to have no particular pattern
anova(model_C_Org_D1)

post_hoc_C_org_D1<-glht(model_C_Org_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_C_org_D1)
#tuk.cld1.C_org <- cld(post_hoc_C_org_D1)

# Date 2
model_C_Org_D2 = lmer(Carbon_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_C_Org_D2),residuals(model_C_Org_D2)) # Needs to have no particular pattern
anova(model_C_Org_D2)

post_hoc_C_org_D2<-glht(model_C_Org_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_C_org_D2)
#tuk.cld2.C_org <- cld(post_hoc_C_org_D2)

# Date 3
model_C_Org_D3 = lmer(Carbon_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_C_Org_D3),residuals(model_C_Org_D3)) # Needs to have no particular pattern
anova(model_C_Org_D3)

post_hoc_C_org_D3<-glht(model_C_Org_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_C_org_D3)
#tuk.cld3.C_org <- cld(post_hoc_C_org_D3)

mean_Corg <-Soil_chem_Org%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(Carbon_PC, na.rm=TRUE))

###### C Organic Boxplot####
Burning_C_evolution_Org <- ggplot(Soil_chem_Org, aes(x=Date,y=Carbon_PC))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(0,35)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Carbon_PC,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Carbon_PC,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  ylab("C (%)")+
  guides(fill = "none")


###### C Mineral MM ####
model_Cmin= lmer(Carbon_PC ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Min,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Cmin),residuals(model_Cmin))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Cmin),residuals(model_Cmin))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Cmin))
# Histogram has to be bell shaped
qqnorm(residuals(model_Cmin))
# Straight line = normal distribution

anova(model_Cmin)

# Date 1
model_C_Min_D1 = lmer(Carbon_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D1,REML=TRUE)
plot(fitted(model_C_Min_D1),residuals(model_C_Min_D1)) # Needs to have no particular pattern
anova(model_C_Min_D1)

post_hoc_C_Min_D1<-glht(model_C_Min_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_C_Min_D1)
#tuk.cld1.C_Min <- cld(post_hoc_C_Min_D1)

# Date 2
model_C_Min_D2 = lmer(Carbon_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D2,REML=TRUE)
plot(fitted(model_C_Min_D2),residuals(model_C_Min_D2)) # Needs to have no particular pattern
anova(model_C_Min_D2)

post_hoc_C_Min_D2<-glht(model_C_Min_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_C_Min_D2)
#tuk.cld2.C_Min <- cld(post_hoc_C_Min_D2)

# Date 3
model_C_Min_D3 = lmer(Carbon_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D3,REML=TRUE)
plot(fitted(model_C_Min_D3),residuals(model_C_Min_D3)) # Needs to have no particular pattern
anova(model_C_Min_D3)

post_hoc_C_Min_D3<-glht(model_C_Min_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_C_Min_D3)
#tuk.cld3.C_Min <- cld(post_hoc_C_Min_D3)


mean_Cmin <-Soil_chem_Min%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(Carbon_PC, na.rm=TRUE))

###### C Mineral Boxplot####

Burning_C_evolution_Min <- ggplot(Soil_chem_Min, aes(x=Date,y=Carbon_PC,color=Burning_intensity))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(0,35)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Carbon_PC,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Carbon_PC,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  ylab("C (%)")+
  guides(fill = "none")

#### Common plot ####

Burning_C_evolution <- ggarrange(Burning_C_evolution_Org,
                                 Burning_C_evolution_Min,
                                 labels= c("A","B"),
                                 ncol=1, 
                                 nrow=2, 
                                 common.legend = TRUE, 
                                 legend="bottom")

#####_______#####

###### N Organic MM ####
hist(Soil_chem_Org$Nitrogen_PC)
 
model_Norg= lmer(Nitrogen_PC ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Org,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Norg),residuals(model_Norg))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Norg),residuals(model_Norg))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Norg))
# Histogram has to be bell shaped
qqnorm(residuals(model_Norg))
# Straight line = normal distribution


anova(model_Norg)

# Date 1
model_N_org_D1 = lmer(Nitrogen_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_N_org_D1),residuals(model_N_org_D1)) # Needs to have no particular pattern
anova(model_N_org_D1)

post_hoc_N_org_D1<-glht(model_N_org_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_N_org_D1)
#tuk.cld1.N_org <- cld(post_hoc_N_org_D1)

# Date 2
model_N_org_D2 = lmer(Nitrogen_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_N_org_D2),residuals(model_N_org_D2)) # Needs to have no particular pattern
anova(model_N_org_D2)

post_hoc_N_org_D2<-glht(model_N_org_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_N_org_D2)
#tuk.cld2.N_org <- cld(post_hoc_N_org_D2)

# Date 3
model_N_org_D3 = lmer(Nitrogen_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_N_org_D3),residuals(model_N_org_D3)) # Needs to have no particular pattern
anova(model_N_org_D3)

post_hoc_N_org_D3<-glht(model_N_org_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_N_org_D3)
#tuk.cld3.N_org <- cld(post_hoc_N_org_D3)

mean_Norg <-Soil_chem_Org%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(Nitrogen_PC, na.rm=TRUE))

###### N  Organic Boxplot ####
Burning_N_evolution_Org  <- ggplot(Soil_chem_Org, aes(x=Date,y=Nitrogen_PC))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(0,1.2)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Nitrogen_PC,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Nitrogen_PC,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  ylab("N (%)")+
  guides(fill = "none")

Burning_N_evolution_Org

###### N Mineral MM ####

hist(Soil_chem_Min$Nitrogen_PC)

model_Nmin= lmer(Nitrogen_PC ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Min,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Nmin),residuals(model_Nmin))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Nmin),residuals(model_Nmin))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Nmin))
# Histogram has to be bell shaped
qqnorm(residuals(model_Nmin))
# Straight line = normal distribution

anova(model_Nmin)

# Date 1
model_N_Min_D1 = lmer(Nitrogen_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D1,REML=TRUE)
plot(fitted(model_N_Min_D1),residuals(model_N_Min_D1)) # Needs to have no particular pattern
anova(model_N_Min_D1)

post_hoc_N_Min_D1<-glht(model_N_Min_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_N_Min_D1)
#tuk.cld1.N_Min <- cld(post_hoc_N_Min_D1)

# Date 2
model_N_Min_D2 = lmer(Nitrogen_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D2,REML=TRUE)
plot(fitted(model_N_Min_D2),residuals(model_N_Min_D2)) # Needs to have no particular pattern
anova(model_N_Min_D2)

post_hoc_N_Min_D2<-glht(model_N_Min_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_N_Min_D2)
#tuk.cld2.N_Min <- cld(post_hoc_N_Min_D2)

# Date 3
model_N_Min_D3 = lmer(Nitrogen_PC ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D3,REML=TRUE)
plot(fitted(model_N_Min_D3),residuals(model_N_Min_D3)) # Needs to have no particular pattern
anova(model_N_Min_D3)

post_hoc_N_Min_D3<-glht(model_N_Min_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_N_Min_D3)
#tuk.cld3.N_Min <- cld(post_hoc_N_Min_D3)


mean_Nmin <-Soil_chem_Min%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(Nitrogen_PC, na.rm=TRUE))


###### N Mineral Boxplot ####
Burning_N_evolution_Min <- ggplot(Soil_chem_Min, aes(x=Date,y=Nitrogen_PC))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(0,1.2)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Nitrogen_PC,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Nitrogen_PC,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  ylab("N (%)")+
  guides(fill = "none")

Burning_N_evolution_Min

#### Common plot ####
Burning_N_evolution <- ggarrange(Burning_N_evolution_Org,
                                 Burning_N_evolution_Min,
                                 labels= c("A","B"),
                                 ncol=1, 
                                 nrow=2, 
                                 common.legend = TRUE, 
                                 legend="bottom")


#####_______#####

###### P Organic MM ####

hist(Soil_chem_Org$P_ppm)#

model_Porg= lmer(P_ppm ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Org,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Porg),residuals(model_Porg))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Porg),residuals(model_Porg))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Porg))
# Histogram has to be bell shaped
qqnorm(residuals(model_Porg))
# Straight line = normal distribution

anova(model_Porg)

mean_Porg<-Soil_chem_Org%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(P_ppm, na.rm=TRUE))

# Date 1 
model_POrg_D1 = lmer(P_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_POrg_D1),residuals(model_POrg_D1)) # Needs to have no particular pattern
anova(model_POrg_D1)

post_hoc_Porg_D1<-glht(model_POrg_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_Porg_D1)
tuk.cld1.Porg <- cld(post_hoc_Porg_D1)


# Date 2 
model_POrg_D2 = lmer(P_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_POrg_D2),residuals(model_POrg_D2)) # Needs to have no particular pattern
anova(model_POrg_D2)

post_hoc_Porg_D2<-glht(model_POrg_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_Porg_D2)
tuk.cld2.Porg <- cld(post_hoc_Porg_D2)


# Date 3 
model_POrg_D3 = lmer(P_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_POrg_D3),residuals(model_POrg_D3)) # Needs to have no particular pattern
anova(model_POrg_D3)

post_hoc_Porg_D3<-glht(model_POrg_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_Porg_D3)
tuk.cld3.Porg <- cld(post_hoc_Porg_D3)



Summarized_Porg <-Soil_chem_Org %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(P_ppm, na.rm=TRUE), Sd = sd(P_ppm,na.rm=T),Max_Porg=max(P_ppm))

dat_text <- data.frame(label = c(tuk.cld1.Porg$mcletters$Letters,tuk.cld2.Porg$mcletters$Letters, tuk.cld3.Porg$mcletters$Letters),
                       Date   = c(rep("2018/06/15",4),rep("2018/09/15",4), rep("2019/09/24",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))


Summarized_Porg_tuk <- merge(Summarized_Porg,dat_text,by=c("Date","Burning_intensity"))

#### P Organic Boxplot ####

Date.labs <- c("Right before burn", "1 month after burn","4 months after burn","16 months after burn")
names(Date.labs) <- c("2018/05/15","2018/06/15", "2018/09/15","2019/09/24")
#Date.labs <- c("1 month after burn","4 months after burn","16 months after burn")
#names(Date.labs) <- c("2018/06/15", "2018/09/15","2019/09/24")

Burning_P_evolution_Org  <- ggplot(Soil_chem_Org, aes(x=Date,y=P_ppm))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(0, 70)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=P_ppm,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=P_ppm,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_Porg_tuk,aes(x=Date,y=2+Max_Porg,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  labs(y=expression(paste("P (",mg.kg^-1,")")))+
  guides(fill = "none")




###### P Mineral MM ####

hist(Soil_chem_Min$P_ppm)

model_Pmin= lmer(P_ppm ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Min,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Pmin),residuals(model_Pmin))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Pmin),residuals(model_Pmin))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Pmin))
# Histogram has to be bell shaped
qqnorm(residuals(model_Pmin))
# Straight line = normal distribution

anova(model_Pmin)

mean_Pmin<-Soil_chem_Min%>%
  group_by(Burning_intensity,Date) %>%
  summarize(Mean = mean(P_ppm, na.rm=TRUE))

# Date 1 
model_PMin_D1 = lmer(P_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D1,REML=TRUE)
plot(fitted(model_PMin_D1),residuals(model_PMin_D1)) # Needs to have no particular pattern
anova(model_PMin_D1)

post_hoc_PMin_D1<-glht(model_PMin_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_PMin_D1)
tuk.cld1.PMin <- cld(post_hoc_PMin_D1)


# Date 2 
model_PMin_D2 = lmer(P_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D2,REML=TRUE)
plot(fitted(model_PMin_D2),residuals(model_PMin_D2)) # Needs to have no particular pattern
anova(model_PMin_D2)


post_hoc_PMin_D2<-glht(model_PMin_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_PMin_D2)
tuk.cld2.PMin <- cld(post_hoc_PMin_D2)

# Date 3 
model_PMin_D3 = lmer(P_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D3,REML=TRUE)
plot(fitted(model_PMin_D3),residuals(model_PMin_D3)) # Needs to have no particular pattern
anova(model_PMin_D3)

post_hoc_PMin_D3<-glht(model_PMin_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_PMin_D3)
tuk.cld3.PMin <- cld(post_hoc_PMin_D3)

Summarized_PMin <-Soil_chem_Min %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(P_ppm, na.rm=TRUE), Sd = sd(P_ppm,na.rm=T),Max_PMin=max(P_ppm))

dat_text <- data.frame(label = c(tuk.cld1.Porg$mcletters$Letters,tuk.cld2.Porg$mcletters$Letters, tuk.cld3.Porg$mcletters$Letters),
                       Date   = c(rep("2018/06/15",4),rep("2018/09/15",4), rep("2019/09/24",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))

Summarized_PMin_tuk <- merge(Summarized_PMin,dat_text,by=c("Date","Burning_intensity"))

###### P  Mineral Boxplot ####

Burning_P_evolution_Min  <- ggplot(Soil_chem_Min, aes(x=Date,y=P_ppm))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(0, 70)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=P_ppm,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=P_ppm,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_PMin_tuk,aes(x=Date,y=1+Max_PMin,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  labs(y=expression(paste("P (",mg.kg^-1,")")))+
  guides(fill = "none")





#### Common plot ####
Burning_P_evolution <- ggarrange(Burning_P_evolution_Org,
                                 Burning_P_evolution_Min,
                                 labels= c("A","B"),
                                 ncol=1, 
                                 nrow=2, 
                                 common.legend = TRUE, 
                                 legend="bottom")



#####_______#####

###### K Organic MM ####
hist(Soil_chem_Org$K_ppm)

model_Korg= lmer(K_ppm ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Org,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Korg),residuals(model_Korg))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Korg),residuals(model_Korg))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Korg))
# Histogram has to be bell shaped
qqnorm(residuals(model_Korg))
# Straight line = normal distribution

anova(model_Korg)

# Date 1
hist(Soil_chem_Org_D1$K_ppm)
model_KOrg_D1 = lmer(K_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_KOrg_D1),residuals(model_KOrg_D1)) # Needs to have no particular pattern
anova(model_KOrg_D1)

post_hoc_KOrg_D1<-glht(model_KOrg_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_KOrg_D1)
tuk.cld1.KOrg <- cld(post_hoc_KOrg_D1)


# Date 2 
model_KOrg_D2 = lmer(K_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_KOrg_D2),residuals(model_KOrg_D2)) # Needs to have no particular pattern
anova(model_KOrg_D2)

post_hoc_KOrg_D2<-glht(model_KOrg_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_KOrg_D2)
tuk.cld2.KOrg <- cld(post_hoc_KOrg_D2)


# Date 3 
model_KOrg_D3 = lmer(K_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_KOrg_D3),residuals(model_KOrg_D3)) # Needs to have no particular pattern
anova(model_KOrg_D3)

post_hoc_KOrg_D3<-glht(model_KOrg_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_KOrg_D3)
tuk.cld3.KOrg <- cld(post_hoc_KOrg_D3)


Summarized_KOrg <-Soil_chem_Org %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(K_ppm, na.rm=TRUE), Sd = sd(K_ppm,na.rm=T),Max_KOrg=max(K_ppm))


dat_text <- data.frame(label = c(tuk.cld1.KOrg$mcletters$Letters,tuk.cld2.KOrg$mcletters$Letters,  tuk.cld3.KOrg$mcletters$Letters),
                       Date   = c(rep("2018/06/15",4),rep("2018/09/15",4), rep("2019/09/24",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))



Summarized_KOrg_tuk <- merge(Summarized_KOrg,dat_text,by=c("Date","Burning_intensity"))

#### K Organic Boxplot ####

Burning_K_evolution_Org  <- ggplot(Soil_chem_Org, aes(x=Date,y=K_ppm))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(0,450)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=K_ppm,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=K_ppm,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_KOrg_tuk,aes(x=Date,y=10+Max_KOrg,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  labs(y=expression(paste("K (",mg.kg^-1,")")))+
  guides(fill = "none")

###### K Mineral MM ####
hist(Soil_chem_Min$K_ppm)
hist(log(Soil_chem_Min$K_ppm))
model_KMin= lmer(log(K_ppm) ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Min,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_KMin),residuals(model_KMin))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_KMin),residuals(model_KMin))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_KMin))
# Histogram has to be bell shaped
qqnorm(residuals(model_KMin))
# Straight line = normal distribution

anova(model_KMin) # Burning effect

# Date 1
hist(Soil_chem_Min_D1$K_ppm)
model_KMin_D1 = lmer(K_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D1,REML=TRUE)
plot(fitted(model_KMin_D1),residuals(model_KMin_D1)) # Needs to have no particular pattern
anova(model_KMin_D1)

post_hoc_KMin_D1<-glht(model_KMin_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_KMin_D1)
tuk.cld1.KMin <- cld(post_hoc_KMin_D1)


# Date 2 
model_KMin_D2 = lmer(K_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D2,REML=TRUE)
plot(fitted(model_KMin_D2),residuals(model_KMin_D2)) # Needs to have no particular pattern
anova(model_KMin_D2)

post_hoc_KMin_D2<-glht(model_KMin_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_KMin_D2)
tuk.cld2.KMin <- cld(post_hoc_KMin_D2)


# Date 3 
model_KMin_D3 = lmer(K_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Min_D3,REML=TRUE)
plot(fitted(model_KMin_D3),residuals(model_KMin_D3)) # Needs to have no particular pattern
anova(model_KMin_D3)

post_hoc_KMin_D3<-glht(model_KMin_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_KMin_D3)
tuk.cld3.KMin <- cld(post_hoc_KMin_D3)

Summarized_KMin <-Soil_chem_Min %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(K_ppm, na.rm=TRUE), Sd = sd(K_ppm,na.rm=T),Max_KMin=max(K_ppm))


dat_text <- data.frame(label = c(tuk.cld1.KMin$mcletters$Letters,tuk.cld2.KMin$mcletters$Letters,  tuk.cld3.KMin$mcletters$Letters),
                       Date   = c(rep("2018/06/15",4),rep("2018/09/15",4), rep("2019/09/24",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))


Summarized_KMin_tuk <- merge(Summarized_KMin,dat_text,by=c("Date","Burning_intensity"))

#### K Mineral Boxplot ####

Burning_K_evolution_Min  <- ggplot(Soil_chem_Min, aes(x=Date,y=K_ppm))+
  geom_boxplot(aes(color=Burning_intensity))+
  ylim(0,450)+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=K_ppm,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=K_ppm,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_KMin_tuk,aes(x=Date,y=10+Max_KMin,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  labs(y=expression(paste("K (",mg.kg^-1,")")))+
  guides(fill = "none")


#### Common plot ####
Burning_K_evolution <- ggarrange(Burning_K_evolution_Org,
                                 Burning_K_evolution_Min,
                                 labels= c("A","B"),
                                 ncol=1, 
                                 nrow=2, 
                                 common.legend = TRUE, 
                                 legend="bottom")

#####_______#####

###### Mg Organic MM ####
hist(Soil_chem_Org$Mg_ppm)
     
model_Mgorg= lmer(Mg_ppm ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Org,REML=TRUE)
     
# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Mgorg),residuals(model_Mgorg))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Mgorg),residuals(model_Mgorg))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Mgorg))
# Histogram has to be bell shaped
qqnorm(residuals(model_Mgorg))
# Straight line = normal distribution

anova(model_Mgorg)

# Date 1
model_MgOrg_D1 = lmer(Mg_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_MgOrg_D1),residuals(model_MgOrg_D1)) # Needs to have no particular pattern
anova(model_MgOrg_D1)

post_hoc_MgOrg_D1<-glht(model_MgOrg_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_MgOrg_D1)
tuk.cld1.MgOrg <- cld(post_hoc_MgOrg_D1)


# Date 2 
model_MgOrg_D2 = lmer(Mg_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_MgOrg_D2),residuals(model_MgOrg_D2)) # Needs to have no particular pattern
anova(model_MgOrg_D2)

post_hoc_MgOrg_D2<-glht(model_MgOrg_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_MgOrg_D2)
tuk.cld2.MgOrg <- cld(post_hoc_MgOrg_D2)


# Date 3 
model_MgOrg_D3 = lmer(Mg_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_MgOrg_D3),residuals(model_MgOrg_D3)) # Needs to have no particular pattern
anova(model_MgOrg_D3)

post_hoc_MgOrg_D3<-glht(model_MgOrg_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_MgOrg_D3)
tuk.cld3.MgOrg <- cld(post_hoc_MgOrg_D3)


Summarized_MgOrg <-Soil_chem_Org %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(Mg_ppm, na.rm=TRUE), Sd = sd(Mg_ppm,na.rm=T),Max_MgOrg=max(Mg_ppm))


dat_text <- data.frame(label = c(tuk.cld1.MgOrg$mcletters$Letters,tuk.cld2.MgOrg$mcletters$Letters,  tuk.cld3.MgOrg$mcletters$Letters),
                       Date   = c(rep("2018/06/15",4),rep("2018/09/15",4), rep("2019/09/24",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))

Summarized_MgOrg_tuk <- merge(Summarized_MgOrg,dat_text,by=c("Date","Burning_intensity"))

#### Mg Organic Boxplot ####

Burning_Mg_evolution_Org  <- ggplot(Soil_chem_Org, aes(x=Date,y=Mg_ppm))+
  ylim(0,450)+
  geom_boxplot(aes(color=Burning_intensity))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Mg_ppm,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Mg_ppm,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_MgOrg_tuk,aes(x=Date,y=25+Max_MgOrg,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  labs(y=expression(paste("Mg (",mg.kg^-1,")")))+
  guides(fill = "none")


###### Mg Mineral MM ####
hist(Soil_chem_Min$Mg_ppm)

model_MgMin <- lmer(Mg_ppm ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Min,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_MgMin),residuals(model_MgMin))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_MgMin),residuals(model_MgMin))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_MgMin))
# Histogram has to be bell shaped
qqnorm(residuals(model_MgMin))
# Straight line = normal distribution

anova(model_MgMin)


# Date 1
model_MgMin_D1 = lmer(Mg_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_MgMin_D1),residuals(model_MgMin_D1)) # Needs to have no particular pattern
anova(model_MgMin_D1)

post_hoc_MgMin_D1<-glht(model_MgMin_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_MgMin_D1)
tuk.cld1.MgMin <- cld(post_hoc_MgMin_D1)


# Date 2 
model_MgMin_D2 = lmer(Mg_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_MgMin_D2),residuals(model_MgMin_D2)) # Needs to have no particular pattern
anova(model_MgMin_D2)

post_hoc_MgMin_D2<-glht(model_MgMin_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_MgMin_D2)
tuk.cld2.MgMin <- cld(post_hoc_MgMin_D2)


# Date 3 
model_MgMin_D3 = lmer(Mg_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_MgMin_D3),residuals(model_MgMin_D3)) # Needs to have no particular pattern
anova(model_MgMin_D3)

post_hoc_MgMin_D3<-glht(model_MgMin_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_MgMin_D3)
tuk.cld3.MgMin <- cld(post_hoc_MgMin_D3)

Summarized_MgMin <-Soil_chem_Min %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(Mg_ppm, na.rm=TRUE), Sd = sd(Mg_ppm,na.rm=T),Max_MgMin=max(Mg_ppm))


dat_text <- data.frame(label = c(tuk.cld1.MgMin$mcletters$Letters,tuk.cld2.MgMin$mcletters$Letters,  tuk.cld3.MgMin$mcletters$Letters),
                       Date   = c(rep("2018/06/15",4),rep("2018/09/15",4), rep("2019/09/24",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))


Summarized_MgMin_tuk <- merge(Summarized_MgMin,dat_text,by=c("Date","Burning_intensity"))

#### Mg Mineral Boxplot ####

Burning_Mg_evolution_Min  <- ggplot(Soil_chem_Min, aes(x=Date,y=Mg_ppm))+
  ylim(0,450)+
  geom_boxplot(aes(color=Burning_intensity))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Mg_ppm,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Mg_ppm,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_MgMin_tuk,aes(x=Date,y=25+Max_MgMin,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  labs(y=expression(paste("Mg (",mg.kg^-1,")")))+
  guides(fill = "none")




#### Common plot ####
Burning_Mg_evolution <- ggarrange(Burning_Mg_evolution_Org,
                                  Burning_Mg_evolution_Min,
                                  labels= c("A","B"),
                                  ncol=1, 
                                  nrow=2, 
                                  common.legend = TRUE, 
                                  legend="bottom")



#####_______#####

###### Ca Organic MM ####
hist(Soil_chem_Org$Ca_ppm)

model_Caorg= lmer(Ca_ppm ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Org,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_Caorg),residuals(model_Caorg))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_Caorg),residuals(model_Caorg))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_Caorg))
# Histogram has to be bell shaped
qqnorm(residuals(model_Caorg))
# Straight line = normal distribution

anova(model_Caorg)

# Date 1
model_CaOrg_D1 = lmer(Ca_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_CaOrg_D1),residuals(model_CaOrg_D1)) # Needs to have no particular pattern
anova(model_CaOrg_D1)

post_hoc_CaOrg_D1<-glht(model_CaOrg_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_CaOrg_D1)
tuk.cld1.CaOrg <- cld(post_hoc_CaOrg_D1)


# Date 2 
model_CaOrg_D2 = lmer(Ca_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_CaOrg_D2),residuals(model_CaOrg_D2)) # Needs to have no particular pattern
anova(model_CaOrg_D2)

post_hoc_CaOrg_D2<-glht(model_CaOrg_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_CaOrg_D2)
tuk.cld2.CaOrg <- cld(post_hoc_CaOrg_D2)


# Date 3 
model_CaOrg_D3 = lmer(Ca_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_CaOrg_D3),residuals(model_CaOrg_D3)) # Needs to have no particular pattern
anova(model_CaOrg_D3)

post_hoc_CaOrg_D3<-glht(model_CaOrg_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_CaOrg_D3)
tuk.cld3.CaOrg <- cld(post_hoc_CaOrg_D3)


Summarized_CaOrg <-Soil_chem_Org %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(Ca_ppm, na.rm=TRUE), Sd = sd(Ca_ppm,na.rm=T),Max_CaOrg=max(Ca_ppm))


dat_text <- data.frame(label = c(tuk.cld1.CaOrg$mcletters$Letters,tuk.cld2.CaOrg$mcletters$Letters,  tuk.cld3.CaOrg$mcletters$Letters),
                       Date   = c(rep("2018/06/15",4),rep("2018/09/15",4), rep("2019/09/24",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))

Summarized_CaOrg_tuk <- merge(Summarized_CaOrg,dat_text,by=c("Date","Burning_intensity"))

#### Ca Organic Boxplot ####

Burning_Ca_evolution_Org  <- ggplot(Soil_chem_Org, aes(x=Date,y=Ca_ppm))+
  ylim(0,3000)+
  geom_boxplot(aes(color=Burning_intensity))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Ca_ppm,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Ca_ppm,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_CaOrg_tuk,aes(x=Date,y=100+Max_CaOrg,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  labs(y=expression(paste("Ca (",mg.kg^-1,")")))+
  guides(fill = "none")


###### Ca Mineral MM ####
hist(Soil_chem_Min$Ca_ppm)

model_CaMin <- lmer(Ca_ppm ~ Burning_intensity * Date + (1|Block/Plot), data=Soil_chem_Min,REML=TRUE)

# Assumptions #
# I. Linearity assumption #
plot(fitted(model_CaMin),residuals(model_CaMin))
# Check if there's a curve / non linear pattern.
# Do a log transformation if there is a pattern

# II. Absence of colineartiy 
par(mfrow=c(1,3))
# III. Homoskedasticity
plot(fitted(model_CaMin),residuals(model_CaMin))
# Needs to have no particular pattern

# IV. Normality of residuals 
hist(residuals(model_CaMin))
# Histogram has to be bell shaped
qqnorm(residuals(model_CaMin))
# Straight line = normal distribution

anova(model_CaMin)


# Date 1
model_CaMin_D1 = lmer(Ca_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D1,REML=TRUE)
plot(fitted(model_CaMin_D1),residuals(model_CaMin_D1)) # Needs to have no particular pattern
anova(model_CaMin_D1)

post_hoc_CaMin_D1<-glht(model_CaMin_D1, mcp(Burning_intensity="Tukey"))
summary(post_hoc_CaMin_D1)
tuk.cld1.CaMin <- cld(post_hoc_CaMin_D1)


# Date 2 
model_CaMin_D2 = lmer(Ca_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D2,REML=TRUE)
plot(fitted(model_CaMin_D2),residuals(model_CaMin_D2)) # Needs to have no particular pattern
anova(model_CaMin_D2)

post_hoc_CaMin_D2<-glht(model_CaMin_D2, mcp(Burning_intensity="Tukey"))
summary(post_hoc_CaMin_D2)
tuk.cld2.CaMin <- cld(post_hoc_CaMin_D2)


# Date 3 
model_CaMin_D3 = lmer(Ca_ppm ~ Burning_intensity + (1|Block), data=Soil_chem_Org_D3,REML=TRUE)
plot(fitted(model_CaMin_D3),residuals(model_CaMin_D3)) # Needs to have no particular pattern
anova(model_CaMin_D3)

post_hoc_CaMin_D3<-glht(model_CaMin_D3, mcp(Burning_intensity="Tukey"))
summary(post_hoc_CaMin_D3)
tuk.cld3.CaMin <- cld(post_hoc_CaMin_D3)



Summarized_CaMin <-Soil_chem_Min %>%
  group_by(Date,Burning_intensity) %>%
  summarize(Mean = mean(Ca_ppm, na.rm=TRUE), Sd = sd(Ca_ppm,na.rm=T),Max_CaMin=max(Ca_ppm))


dat_text <- data.frame(label = c(tuk.cld1.CaMin$mcletters$Letters,tuk.cld2.CaMin$mcletters$Letters,  tuk.cld3.CaMin$mcletters$Letters),
                       Date   = c(rep("2018/06/15",4),rep("2018/09/15",4), rep("2019/09/24",4)),
                       Burning_intensity   = rep(c("NoBurn","LowBurn","MedBurn","HighBurn"),3))

Summarized_CaMin_tuk <- merge(Summarized_CaMin,dat_text,by=c("Date","Burning_intensity"))

#### Ca Mineral Boxplot ####

Burning_Ca_evolution_Min  <- ggplot(Soil_chem_Min, aes(x=Date,y=Ca_ppm))+
  ylim(0,3000)+
  geom_boxplot(aes(color=Burning_intensity))+
  facet_grid(.~Date,scales="free_x",drop=T, labeller=labeller(Date=Date.labs))+
  stat_summary(fun="mean", geom="point",shape=21,size=2,aes(x=Date,y=Ca_ppm,fill=Burning_intensity),
               position=position_dodge(width=0.75))+ # Add a black contour over mean, position dodge 0.75 matches the dodge with boxplot
  stat_summary(fun="mean", geom="point",size=1.3,aes(x=Date,y=Ca_ppm,color=Burning_intensity),
               position=position_dodge(width=0.75))+
  geom_text(data=Summarized_CaMin_tuk,aes(x=Date,y=200+Max_CaMin,label=label,group=Burning_intensity),position=position_dodge(width=0.75))+
  theme_light()+
  theme(text =  element_text(size=15))+
  theme(legend.position = "bottom",
        legend.background = element_rect(size=0.4, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=16,face="bold"))+
  scale_colour_manual(name= "Burning intensity",values = c("#9cccd8","#ffb53f","#f26419","#931125"))+
  theme(strip.text.x = element_text(size = 12, colour = "black",face="bold"))+
  theme(axis.title.x = element_blank())+
  labs(y=expression(paste("Ca (",mg.kg^-1,")")))+
  guides(fill = "none")




#### Common plot ####
Burning_Ca_evolution <- ggarrange(Burning_Ca_evolution_Org,
                                  Burning_Ca_evolution_Min,
                                  labels= c("A","B"),
                                  ncol=1, 
                                  nrow=2, 
                                  common.legend = TRUE, 
                                  legend="bottom")
Burning_Ca_evolution




#####_______#####
Burning_pH_evolution
Burning_C_evolution
Burning_N_evolution
Burning_P_evolution
Burning_K_evolution
Burning_Ca_evolution
Burning_S_evolution
Burning_Mg_evolution
