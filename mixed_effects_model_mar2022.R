# Corinne Klohmann
# cak268@uw.edu

# Mixed effects GLM with negative binomial distribution 

# clear workspace
rm(list =ls())

## load packages
require(foreign)
require(ggplot2)
require(MASS)
require(lme4)
# load data
dat <- read.csv("Field Assay Counts - Sheet1.csv")

# clean up data
#remove discovery field site
dat <- dat[dat$site != "Disc",]
# remove geom mean, notes columns
dat <- dat[-c(15,16, 17,18)]
# remove sample
dat <- dat[,-12]

#see data
summary(dat)
# make date julian day
library(lubridate)
# convert "date" from chr to a Date class and specify current date format
dat$Date<- as.Date(dat$Date, "%m/%d/%y")
# convert with yday into a new column "julian"
dat$julian <- yday(dat$Date)  

# make sure it worked all the way through. 
head(dat$julian) 


tail(dat$julian)

# model
m.glm <- glm.nb(ent_count ~ grass.nograss + dist.from.shore + do,
                data=dat, trace=TRUE)
summary(m.glm)

## compare models 

ent.nb1 <- glm.nb(ent_count ~ grass.nograss + dist.from.shore + do,
                  data=dat, trace = TRUE) ## best model lowest AIC
ent.nb2 <- update(ent.nb1, ent_count ~ grass.nograss + dist.from.shore + do + dist.from.shore*do)
ent.nb3 <- update(ent.nb2, ent_count ~ grass.nograss + dist.from.shore)
anova(ent.nb1, ent.nb2, ent.nb3)
# rescale predictor variables 

dat$dist.from.shore<-scale(dat$dist.from.shore)
dat$do<-scale(dat$do)
dat$grass.nograss<- scale(dat$grass.nograss)

# add in random effects (date and site)
m.nb <- glmer.nb(ent_count ~ grass.nograss + dist.from.shore + do + (1|julian) + (1|site),
                 data=dat, verbose=TRUE)
m.nb

#plotting predicted values
plot(m.nb)
summary(m.nb)



# remove nas
dat <- na.omit(dat)


  
##### visualize results 
summ(m.nb)

# MODEL INFO:
#   Observations: 712
# Dependent Variable: ent_count
# Type: Mixed effects generalized linear regression
# Error Distribution: Negative Binomial(1.5516)
# Link function: log 
# 
# MODEL FIT:
#   AIC = 5961.98, BIC = 5993.96
# Pseudo-R² (fixed effects) = 0.03
# Pseudo-R² (total) = 0.97 
# 
# FIXED EFFECTS:
#   ----------------------------------------------------
#   Est.   S.E.   z val.      p
# --------------------- ------- ------ -------- ------
#   (Intercept)              3.11   0.47     6.56   0.00
# grass.nograss           -0.00   0.04    -0.07   0.95
# dist.from.shore         -0.16   0.03    -4.67   0.00
# do                       0.13   0.11     1.13   0.26
# ----------------------------------------------------
#   
#   RANDOM EFFECTS:
#   ----------------------------------
#   Group     Parameter    Std. Dev. 
# -------- ------------- -----------
#   julian   (Intercept)     0.92    
# site    (Intercept)     0.71    
# ----------------------------------
#   
#   Grouping variables:
#   --------------------------
#   Group    # groups   ICC  
# -------- ---------- ------
#   julian      17      0.25 
# site       3       0.12 
# --------------------------
# explore the impact of dist from shore on ent abundance 
effect_plot(m.nb, pred = dist.from.shore)
# uncertainty around model predictions
effect_plot(m.nb, pred = dist.from.shore, interval = TRUE)
# how the data are distributed
effect_plot(m.nb, pred = dist.from.shore, interval = TRUE, rug = TRUE)
# model relates to the observed data
effect_plot(m.nb, pred = dist.from.shore, interval = TRUE, plot.points = TRUE)
effect_plot(m.nb, pred = do, interval = TRUE, plot.points = TRUE)
# model relates to the observed data
effect_plot(m.nb, pred = do, interval = TRUE, plot.points = TRUE)
# Partial residuals plots
effect_plot(m.nb, pred = dist.from.shore, interval = TRUE, partial.residuals = TRUE)
effect_plot(m.nb, pred = do, interval = TRUE, partial.residuals = TRUE)
# add jigger 
effect_plot(m.nb, pred = dist.from.shore, interval = TRUE, partial.residuals = TRUE,
            jitter = c(0.2,0))


# categorical predictors 
effect_plot(m.nb, pred = grass.nograss, interval = TRUE)
effect_plot(m.nb, pred = grass.nograss, interval = TRUE, plot.points = TRUE,
            jitter = .2)
effect_plot(m.nb, pred = grass.nograss, interval = TRUE, partial.residuals = TRUE,
            jitter = .2)

# plot 
library(effects) 
plot(allEffects(m.nb))
plot(allEffects(m.nb), ylab = "Enterococcus (CFUs/100ml)")

plot(m.nb, which=c(1,5))

m.nb2 <- allEffects(m.nb)
plot(m.nb2[4], x.var="dist.from.shore", multiline=TRUE, ci.style="bands", ...)
plot(.mnb2[4], multiline=TRUE, ci.style="bands", ...)

# plot predictors
ent_plot <- ggplot(dat,
                    aes(x = dist.from.shore, y = ent_count )) +
  geom_point()

ent_plot
# plot residuals 
library(tidyverse)
library(readr)
library(readxl)
# make plot - fits well!
m.nb_res <- simulateResiduals(m.nb)
plot(m.nb_res)







##### new try ######  #### GOOD ONE ###  
library(tidyverse) #for all data wrangling
library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions

## Unformatted plot of effect sizes
sjPlot::plot_model(m.nb)
## Formatted plot of effect sizes
# Notes: axis labels should be in order from bottom to top. 
# To see the values of the effect size and p-value, set show.values and show.p= TRUE. Pvalues will only be shown if the effect size values are too

b <-sjPlot::plot_model(m.nb, 
                   axis.labels=c("Dissolved Oxygen (Mg/L)", "Distance from shore (m)", 
                                 "Seagrass"),
                   show.values=TRUE, show.p=TRUE,
                   title= "Field effects on Enterococcus") +
  theme(panel.grid = element_blank())
# unformatted table
sjPlot:: tab_model(m.nb)
# Formatted table
sjPlot::tab_model(m.nb, 
                  show.re.var= TRUE, 
                  pred.labels =c("(Intercept)", "Z. marina", "Distance from shore (m)", 
                                 "Dissolved oxygen (Mg/L)"),
                  dv.labels= "Field effects on Enterococcus Abundance")
## Plot model estimates WITH data ##
# Step 1: Save the effect size estimates into a data.frame
# Use the effects package --> effect function
# term= the fixed effect you want to get data on, mod= name of your model.

effects_m.nb <- effects::effect(term= "dist.from.shore", mod= m.nb)
summary(effects_m.nb) #output of what the values are
# Save the effects values as a df:
x_m.nb <- as.data.frame(effects_m.nb)
x_m.nb$dist.from.shore <- c(0,25,50,75,100)

# Step 2: Use the effects value df (created above) to plot the estimates
#1
m.nb_plot <- ggplot() + 
  #2
  geom_point(data=dat, aes(dist.from.shore, ent_count)) + 
  #3
  geom_point(data=x_m.nb, aes(x=dist.from.shore, y=fit), color="blue") +
  #4
  geom_line(data=x_m.nb, aes(x=dist.from.shore, y=fit), color="blue") +
  #5
  geom_ribbon(data= x_m.nb, aes(x=dist.from.shore, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  #6
  labs(x="Distance from Shore (m)", y="Enterococcus (CFUs/100ml)")

 a <- m.nb_plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 14))
 # packages 
 install.packages("ggpubr")
 library(ggpubr)
#combine plots 
figure <- ggarrange(a, b,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)


# Step 1: Save the effect size estimates into a data.frame
# Use the effects package --> effect function
# term= the fixed effect you want to get data on, mod= name of your model.

effects_m <- effects::effect(term= "do", mod= m.nb)
summary(effects_m) #output of what the values are
# Save the effects values as a df:
x_m <- as.data.frame(effects_m)

# Step 2: Use the effects value df (created above) to plot the estimates
#1
m_plot <- ggplot() + 
  #2
  geom_point(data=dat, aes(do, ent_count)) + 
  #3
  geom_point(data=x_m, aes(x=do, y=fit), color="blue") +
  #4
  geom_line(data=x_m, aes(x=do, y=fit), color="blue") +
  #5
  geom_ribbon(data= x_m, aes(x=do, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  #6
  labs(x="Dissolved Oxygen (Mg/L)", y="Enterococcus (CFUs/100ml)")

m_plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Step 1: Save the effect size estimates into a data.frame
# Use the effects package --> effect function
# term= the fixed effect you want to get data on, mod= name of your model.

effects_e <- effects::effect(term= "grass.nograss", mod= m.nb)
summary(effects_e) #output of what the values are
# Save the effects values as a df:
x_e <- as.data.frame(effects_e)

# Step 2: Use the effects value df (created above) to plot the estimates
#1
m_plot <- ggplot() + 
  #2
  geom_point(data=dat, aes(grass.nograss, ent_count)) + 
  #3
  geom_point(data=x_e, aes(x=grass.nograss, y=fit), color="blue") +
  #4
  geom_line(data=x_e, aes(x=grass.nograss, y=fit), color="blue") +
  #5
  geom_ribbon(data= x_e, aes(x=grass.nograss, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  #6
  labs(x="Grass presence", y="Enterococcus (CFUs/100ml)")

m_plot












