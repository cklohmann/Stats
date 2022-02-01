# Corinne Klohmann 
# Dec. 2021
# cak268@uw.edu

# Generlized linear model for field data 

# clear workspace
rm(list = ls())

# Load required packages
library(googledrive)
library(skimr)
library(lubridate)
library(tidyverse)
library(clm)
library(ordinal)


# Verify that tidyverse is up to date
tidyverse_update()

# read in the data
fdata <- read.csv("Field Assay Counts - Sheet1.csv")

# Count everything up

port_data<-fdata %>%
  filter(ent_count >= 0)
length(port_data$ent_count)

# How many samples in the full dataset?
n_full <- fdata %>%
  group_by(grass.nograss, site) %>%
  summarise(n = n())




########### PORTUNION ########### 

port_model <- ent_count ~ grass.nograss*site + Date 

# Fit a normal distribution
model_norm_port<-glm(formula = port_model, data=fdata)
summary(model_norm_port)	
AIC(model_norm_port)
#AIC = 8355.74
#SELECT THIS MODEL (lowest AIC)

model_pois_port<-glm(formula = port_model, data=fdata, family="poisson")
summary(model_pois_port)	
AIC(model_pois_port)
#AIC = 17441.57

# this model not working, glm.nb() not recognized
model_nbinom_port<-glm.nb(formula = port_model, data=fdata)
summary(model_nbinom_port)	
AIC(model_nbinom_port)
#AIC = 207.1491

expt_port<-ggpredict(model_norm_port,c("grass.nograss", "site"))

## now plot

portunion_plot<-ggplot(expt_port,aes(as.factor(x),predicted,color=group ),
          grouping=group,color=group) + geom_point(size=4) +
 geom_errorbar(data=expt_port,mapping=aes(x=x,ymin=conf.low,ymax=conf.high),width=0.03)+
 geom_line(aes(group=group))+
 xlab("")#
  #ylab(expression(paste("predicted abundance of ",italic("P. conformis"))))+
 #theme_minimal()+
 #theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=9),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
 #scale_x_discrete(limits=rev(levels(expt_port$x)),labels=c("grass","no grass"))+
 #theme(legend.position="none")
portunion_plot
