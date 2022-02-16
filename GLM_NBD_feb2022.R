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
dat <- dat[-c(14,15)]
# remove Date, sample
dat <- dat[-c(1,11)]

#see data
summary(dat)

# plot the data
ggplot(dat, aes(ent_count, fill = site)) + geom_histogram(binwidth = 1) +
  facet_grid(site ~ ., margins = TRUE, scales = "free")

#GLM analysis 
summary(m1 <- glm.nb(ent_count ~ grass.nograss + site, data = dat))
## Call:
# glm.nb(formula = ent_count ~ site + season, data = dat, init.theta = 1.032713156, 
 #   link = log)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -2.5390  -1.0969  -0.437   0.2403   4.7314  
## 
## Coefficients:
##                Estimate Std. Error z value Pr(>|z|)    
## (Intercept)     1.68633    0.09027   18.681  < 2e-16 ***
## grass.nograss   0.17674    0.08397   2.105    0.0353 *  
## siteGG          2.05429    0.10514   19.538   <2e-16 ***  
## siteSmith's Cove  2.19663    0.10616  20.691   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for Negative Binomial(1.033) family taken to be 1)
## 
##     Null deviance: 1331.46  on 806  degrees of freedom
## Residual deviance:  948.16  on 803  degrees of freedom
## (21 observations deleted due to missingness)
## AIC: 6894.7

## Number of Fisher Scoring iterations: 1


## Theta:  0.7362 
## Std. Err.:  0.0366 

## 2 x log-likelihood:  -6884.7330 

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
# 2 -8.094503e-11       1

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



