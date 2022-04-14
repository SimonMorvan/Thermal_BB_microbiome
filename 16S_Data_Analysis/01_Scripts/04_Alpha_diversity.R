# Simon Morvan
# Sept 2O21
R.version # 4.1.1.

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                               Alpha diversity script                       #                                   
#                                                                            #                            #
#                                                                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

Dir <- "/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/16S_Data_Analysis/"
setwd(Dir)

#### Packages ####
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggrepel)

#### Data ####
load("03_RData/ps_bac_merg.RData")

#### Non-rarefied ####

alpha.div<- plot_richness(ps.bac.merg, x="Burning_intensity", measures=c("Simpson", "Shannon"))

# SIMPSON
Simpson <- ggplot(alpha.div$data[alpha.div$data$variable %in% 'Simpson',], aes(x=Burning_intensity, y=value))+
  geom_point(colour = "black", size = 4.5) + 
  geom_point(aes(color=Burning_intensity),size=4)+
  geom_text_repel(aes(label=Plot), color="black")+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),colour="azure4")+
  ggtitle("Simpson reciprocal index")+
  ylab("Alpha diversity index")+
  theme(title = element_text(size=20))+
  theme(axis.text.x =element_blank(),axis.text.y =element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=15))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=20))+
  theme(legend.position = "bottom")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))


# SHANNON 
Shannon <- ggplot(alpha.div$data[alpha.div$data$variable %in% 'Shannon',], aes(x=Burning_intensity, y=value))+
  geom_point(colour = "black", size = 4.5) + 
  geom_point(aes(color=Burning_intensity),size=4)+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),colour="azure4")+
  geom_text_repel(aes(label=Plot), color="black")+
  ggtitle(" Shannon-Weaver Index (H')")+
  ylab("Alpha diversity index")+
  theme(title = element_text(size=20))+
  theme(axis.text.x =element_blank(),axis.text.y =element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=15))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=20))+
  theme(legend.position = "bottom")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))

ggarrange(Shannon, Simpson, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

# Verification with vegan
# Simpson.index <- alpha.div.Plot.asv$data[alpha.div.Plot.asv$data$variable %in% 'Simpson',]
# Shannon.div <- alpha.div.Plot.asv$data[alpha.div.Plot.asv$data$variable %in% 'Shannon',]
# 
# Simps.div.veg <- as.data.frame(diversity(t(ASV.noSD), index="simpson"))
# Shan.div.veg <- as.data.frame(diversity(t(ASV.noSD), index="shannon"))
# The indices found found by phyloseq are the same found by vegan.


####  ANOVA #####
####  SHANNON 
Shannon.div <- alpha.div$data[alpha.div$data$variable %in% 'Shannon',]

res.aov <- aov(value ~ Burning_intensity + Block, data = Shannon.div)
summary(res.aov)

# ANOVA assumptions #
# 1. Variance homogeneity 
# If p-value > 0.05 variance is homogeneous 
bartlett.test(value~Burning_intensity, data=Shannon.div) 

# 2. Residuals normality
#  Q-Q plot 
plot(res.aov, 2)  
# Points have to be on the line

# Shapiro-Wilk tet
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )
# f p-value > 0.05 residuals are normaly distributed.

# 3. Residuals homoscedacity
plot(aov_residuals)
# Needs to have random pattern

# Post-hoc test to see whichgroup have signif. dif. 
Tukey.test <- TukeyHSD(res.aov, ordered = TRUE)


####  SIMPSON 
Simpson.div <- alpha.div$data[alpha.div$data$variable %in% 'Simpson',]

res.aov <- aov(value ~ Burning_intensity + Block, data = Simpson.div)
summary(res.aov)

# ANOVA assumptions #
# 1. Variance homogeneity 
# If p-value > 0.05 variance is homogeneous 
bartlett.test(value~Burning_intensity, data=Simpson.div) 

# 2. Residuals normality
#  Q-Q plot 
plot(res.aov, 2)  

#Shapiro-Wilk test
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

# 3. Residuals homoscedacity
plot(aov_residuals)

Tukey.test <- TukeyHSD(res.aov, ordered = TRUE)



#### Rarefied #####
set.seed(124)
ps.bac.rare <- rarefy_even_depth(ps.bac.merg, sample.size = min(sample_sums(ps.bac.merg)),
                                 rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

# 2 ASV removed.


alpha.div.rar <- plot_richness(ps.bac.rare, x="Burning_intensity", measures=c("Simpson", "Shannon"))

# SIMPSON
Simpson.rar.bac <- ggplot(alpha.div.rar$data[alpha.div.rar$data$variable %in% 'Simpson',], aes(x=Burning_intensity, y=value))+
  geom_point(colour = "black", size = 4.5) + 
  geom_point(aes(color=Burning_intensity),size=4)+
  geom_text_repel(aes(label=Plot), color="black")+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),colour="azure4")+
  ggtitle("Simpson reciprocal index")+
  ylab("Alpha diversity index")+
  theme(title = element_text(size=20))+
  theme(axis.text.x =element_blank(),axis.text.y =element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=15))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=20))+
  theme(legend.position = "bottom")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))

# SHANNON 
Shannon.rar.bac <- ggplot(alpha.div.rar$data[alpha.div.rar$data$variable %in% 'Shannon',], aes(x=Burning_intensity, y=value))+
  geom_point(colour = "black", size = 4.5) + 
  geom_point(aes(color=Burning_intensity),size=4)+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),colour="azure4")+
  geom_text_repel(aes(label=Plot), color="black")+
  ggtitle(" Shannon-Weaver Index (H')")+
  ylab("Alpha diversity index")+
  theme(title = element_text(size=20))+
  theme(axis.text.x =element_blank(),axis.text.y =element_text(size=12))+ # size of lables, axis titles
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=15))+
  theme(axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=20))+
  theme(legend.position = "bottom")+
  theme(legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="black"))+
  scale_color_manual(values = c("#9cccd8","#ffb53f","#f26419","#931125"))


ggarrange(Shannon.rar.bac, Simpson.rar.bac, 
          ncol=2, nrow=1, 
          common.legend = TRUE, legend="bottom")


Shannon.div.rar <- alpha.div.rar$data[alpha.div.rar$data$variable %in% 'Shannon',]
Simpson.div.rar <- alpha.div.rar$data[alpha.div.rar$data$variable %in% 'Simpson',]



# ANOVA 
res.aov <- aov(value ~ Burning_intensity + Block, data = Shannon.div.rar)
summary(res.aov)

# ANOVA assumptions #
# 1. Variance homogeneity 
# If p-value > 0.05 variance is homogeneous 
bartlett.test(value~Burning_intensity, data=Shannon.div.rar) 

# 2. Residuals normality
#  Q-Q plot 
plot(res.aov, 2)  

# Shapiro Wilk test
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

# 3. Residuals homoscedacity
plot(aov_residuals)

Tukey.test <- TukeyHSD(res.aov, ordered = TRUE)


# ANOVA 
res.aov <- aov(value ~ Burning_intensity + Block, data = Simpson.div.rar)
summary(res.aov)

# ANOVA assumptions #
# 1. Variance homogeneity 
# If p-value > 0.05 variance is homogeneous 
bartlett.test(value~Burning_intensity, data=Simpson.div.rar) 

# 2. Residuals normality
#  Q-Q plot 
plot(res.aov, 2) 

# Shapiro-Wilks test
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals )

# 3. Residuals homoscedacity
plot(aov_residuals)

Tukey.test <- TukeyHSD(res.aov, ordered = TRUE)




