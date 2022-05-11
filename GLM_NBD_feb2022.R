# Corinne Klohmann
# cak268@uw.edu

# generalized linear model witn neg biom dist

# clear workspace
rm(list =ls())

# load packages
require(foreign)
require(ggplot2)
require(MASS)

# load data
dat <- read.csv("Field Assay Counts - Sheet1.csv")

# clean up data
#remove discovery field site
dat <- dat[dat$site != "Disc",]
# remove geom mean, notes columns
dat <- dat[-c(15,16, 17,18)]
# remove Date, sample
dat <- dat[-c(1,12)]

#see data
summary(dat)

# correlation table
# corrplot()

# plot the data
ggplot(dat, aes(ent_count, fill = dist.from.shore)) + geom_histogram(binwidth = 1) +
  facet_grid(dist.from.shore ~ ., margins = TRUE, scales = "free")
#GLM analysis 
summary(m1 <- glm.nb(ent_count ~ grass.nograss + dist.from.shore + 
                       do + do*dist.from.shore , data = dat))
# Call:
#   glm.nb(formula = ent_count ~ grass.nograss + dist.from.shore + 
#            do + do * dist.from.shore, data = dat, init.theta = 0.6191299819, 
#          link = log)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3896  -1.1705  -0.4152   0.2855   2.5853  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         5.4110046  0.5120066  10.568  < 2e-16 ***
#   grass.nograss       0.2137385  0.0979042   2.183 0.029026 *  
#   dist.from.shore    -0.0013167  0.0081083  -0.162 0.871002    
# do                 -0.2153352  0.0615273  -3.500 0.000466 ***
#   dist.from.shore:do -0.0001691  0.0009742  -0.174 0.862191    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(0.6191) family taken to be 1)
# 
# Null deviance: 902.33  on 711  degrees of freedom
# Residual deviance: 851.28  on 707  degrees of freedom
# (116 observations deleted due to missingness)
# AIC: 6531.1
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  0.6191 
# Std. Err.:  0.0309 
# 
# 2 x log-likelihood:  -6519.0850 

##### do*dist.from.shore is not significant and the AIC value is not the smallest ##
#### so we drop this term from the model ######


#GLM analysis 
summary(m1 <- glm.nb(ent_count ~ grass.nograss + dist.from.shore + 
                       do, data = dat))
# Call:
#   glm.nb(formula = ent_count ~ grass.nograss + dist.from.shore + 
#            do, data = dat, init.theta = 0.6191012584, link = log)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3920  -1.1690  -0.4186   0.2835   2.6250  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      5.483887   0.282374  19.421  < 2e-16 ***
#   grass.nograss    0.215010   0.097588   2.203   0.0276 *  
#   dist.from.shore -0.002716   0.001429  -1.900   0.0574 .  
# do              -0.224127   0.032752  -6.843 7.75e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(0.6191) family taken to be 1)
# 
# Null deviance: 902.30  on 711  degrees of freedom
# Residual deviance: 851.29  on 708  degrees of freedom
# (116 observations deleted due to missingness)
# AIC: 6529.1
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  0.6191 
# Std. Err.:  0.0309 
# 
# 2 x log-likelihood:  -6519.1240  

#GLM analysis without DO 
summary(m1 <- glm.nb(ent_count ~ grass.nograss + dist.from.shore, data = dat))
# Call:
#   glm.nb(formula = ent_count ~ grass.nograss + dist.from.shore, 
#          data = dat, init.theta = 0.5168934414, link = log)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.1354  -1.2356  -0.4639   0.2494   3.0923  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      3.637647   0.102707  35.418   <2e-16 ***
#   grass.nograss    0.102934   0.098734   1.043    0.297    
# dist.from.shore -0.002003   0.001460  -1.372    0.170    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(0.5169) family taken to be 1)
# 
# Null deviance: 973.57  on 806  degrees of freedom
# Residual deviance: 970.48  on 804  degrees of freedom
# (21 observations deleted due to missingness)
# AIC: 7209.6
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  0.5169 
# Std. Err.:  0.0241 
# 
# 2 x log-likelihood:  -7201.6290 

# the difference in AIC values
6529.1 - 7209.6 # -680.5 (so model with DO is better, no autocorelation!)



# see if site is significant
m2 <- update(m1, . ~ . - site)
anova(m1, m2)


## Likelihood ratio tests of Negative Binomial Models

## Response: ent_count
## Model     theta Resid. df    2 x log-lik.   Test    df
## 1        grass.nograss 0.5158198       805       -7203.575             
## 2 grass.nograss + site 0.7361638       803       -6884.733 1 vs 2     2
## LR stat. Pr(Chi)
## 1                 
## 2 318.8426       0   <- NOT significant predictor or ent

# see if specific sites are significant
# GG 
m2 <- update(m1, . ~ . - siteGG)
anova(m1, m2)

# Likelihood ratio tests of Negative Binomial Models

# Response: ent_count
# Model     theta Resid. df    2 x log-lik.   Test    df
# 1 grass.nograss + site 0.7361638       803       -6884.733             
# 2 grass.nograss + site 0.7361638       803       -6884.733 1 vs 2     0
# LR stat. Pr(Chi)
# 1                      
# 2 -8.094503e-11       1     <- GG IS a significant predictor 

## Alki
m2 <- update(m1, . ~ . - siteAlki)
anova(m1, m2)
## Likelihood ratio tests of Negative Binomial Models

## Response: ent_count
# Model     theta Resid. df    2 x log-lik.   Test    df
# 1 grass.nograss + site 0.7361638       803       -6884.733             
# 2 grass.nograss + site 0.7361638       803       -6884.733 1 vs 2     0
# LR stat. Pr(Chi)
# 1                      
# 2 -8.094503e-11       1   <- Alki IS significant predictor 

## Smith's Cove
m2 <- update(m1, . ~ . - sitesmithcove)
anova(m1, m2)

## Likelihood ratio tests of Negative Binomial Models

## Response: ent_count
## 1 grass.nograss + site 0.7361638       803       -6884.733             
## 2 grass.nograss + site 0.7361638       803       -6884.733 1 vs 2     0
## LR stat. Pr(Chi)
# 1                      
# 2 -8.094503e-11       1  <- SC IS significant predictor 

###### Checking model assumptions ####
# estimating a dispersion parameter (not shown in the output) 
# that is held constant in a Poisson model
m3 <- glm(ent_count ~ grass.nograss + site, family = "poisson", data = dat)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)

## 'log Lik.' 0 (df=5)


## get the confidence intervals for the coefficients 
## by profiling the likelihood function
(est <- cbind(Estimate = coef(m1), confint(m1)))
## Waiting for profiling to be done...
##                Estimate      2.5 %    97.5 %
## (Intercept)   1.6863303 1.50835497 1.8702118
## grass.nograss 0.1767446 0.01153247 0.3417719
## siteGG        2.0542917 1.84754610 2.2600821
## sitesmithcove 2.1966251 1.98749543 2.4050808

# Look at incident rate ratios 
# exponentiate our model coefficients. 
# same applies to the confidence intervals.

exp(est)

##             Estimate    2.5 %    97.5 %
## (Intercept)   5.399629 4.519290  6.489671
## grass.nograss 1.193326 1.011599  1.407439
## siteGG        7.801311 6.344232  9.583876
## sitesmithcove 8.994607 7.297234 11.079326

# predicted values - can't do this with site and grass.nograss
# make site a factor
dat$site <- as.factor(dat$site)
#run
newdata1 <- data.frame(grass.nograss = mean(dat$grass.nograss), 
                       site = factor(1:3, levels = 1:3, 
                    labels = levels(dat$site)))
newdata1$phat <- predict(m1, newdata1, type = "response")
newdata1

## grass.nograss      site      phat
## 1      0.513285      Alki  5.912393
## 2      0.513285        GG 46.124416
## 3      0.513285 smithcove 53.179651

## obtain the mean predicted number of events for values of grass no grass
# across its entire range for each level of site and graph these

newdata2 <- data.frame(
  math = rep(seq(from = min(dat$grass.nograss), to = max(dat$grass.nograss), 
                 length.out = 100), 3),
  prog = factor(rep(1:3, each = 100), levels = 1:3, labels =
                  levels(dat$site)))

newdata2 <- cbind(newdata2, predict(m1, newdata2, type = "link", se.fit=TRUE))
newdata2 <- within(newdata2, {
  DaysAbsent <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata2, aes(math, DaysAbsent)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = prog), alpha = .25) +
  geom_line(aes(colour = prog), size = 2) +
  labs(x = "Grass present or absent", y = "Predicted Ent Count")



