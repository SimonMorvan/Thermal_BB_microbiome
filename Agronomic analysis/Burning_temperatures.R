#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 
#     Visualization of the temperatures recorded at 1cm depth for the 
#                        different burning intensities 
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


##### PACKAGES 
library(reshape2)
library(ggplot2)
library(plyr)
library(directlabels)
library(stringr)
library(ggpubr)

#### DATA ####
setwd("/Users/simonmorvan/Desktop/UdeM/Projet/Ch2/Data/")
Burning_1cm <- read.csv2("Chaleur_1cm.csv")
str(Burning_1cm)
Burning_1cm <-  Burning_1cm[,-1]
Burning_1cm<-lapply(Burning_1cm,as.numeric)
Burning_1cm <- as.data.frame(Burning_1cm)
Depth <- melt(Burning_1cm,id.vars=c("Minutes"))

Depth$intensity <- Depth$variable
revalue(Depth$intensity, c("Plot_3" = "MedBurn","Plot_4" = "LowBurn","Plot_7" = "LowBurn","Plot_8" = "HighBurn","Plot_11" = "MedBurn","Plot_12" = "HighBurn","Plot_15" = "HighBurn","Plot_16" = "MedBurn")) -> Depth$intensity
Depth$intensity<- factor(Depth$intensity , levels = c( "LowBurn", "MedBurn","HighBurn"))


##### FIGURE S1 #####
Depth_1cm <- ggplot(Depth, aes(x=Minutes,y=value,color=intensity))+
             geom_path()+
             geom_dl(aes(label = variable), method = list(dl.combine("first.points", "last.points"),cex = 0.6),color="black") +
             scale_colour_manual(name= "Burning intensity",values = c("#ffb53f","#f26419","#931125"))+
             ylab("Temperature (°C)")+
             xlab("Time since burn (min)")+
             geom_vline(xintercept = 0, color="Black", linetype =2)+
             theme(legend.position = "bottom",
                   legend.background = element_rect(size=0.4, 
                                                    linetype="solid", 
                                                    colour ="black"),
                   legend.text = element_text(size=15),
                   legend.title = element_text(size=16,face="bold"))


