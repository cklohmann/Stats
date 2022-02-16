# Corinne Klohmann
# cak268@uw.edu

# generalized linear model witn neg biom dist
# site and DO

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
# make site a factor
dat$grass.nograss <- as.factor(dat$grass.nograss)
#remove NA's
dat <- na.omit(dat)
#see data
summary(dat)

# plot the data
ggplot(dat, aes(ent_count, fill = site)) + geom_histogram(binwidth = 1) +
  facet_grid(site ~ ., margins = TRUE, scales = "free")

#GLM analysis 
summary(m1 <- glm.nb(ent_count ~ grass.nograss + do, data = dat))
# Call:
#   glm.nb(formula = ent_count ~ grass.nograss + do, data = dat, 
#          init.theta = 0.5878148481, link = log)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3582  -1.1742  -0.4567   0.2816   2.7404  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     5.37816    0.28208  19.066  < 2e-16 ***
#   grass.nograss1  0.26389    0.10291   2.564   0.0103 *  
#   do             -0.23195    0.03518  -6.593 4.32e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(0.5878) family taken to be 1)
# 
# Null deviance: 855.21  on 675  degrees of freedom
# Residual deviance: 809.95  on 673  degrees of freedom
# AIC: 6193
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  0.5878 
# Std. Err.:  0.0300 
# 
# 2 x log-likelihood:  -6184.9970  

# see if DO is significant
m2 <- update(m1, . ~ . - do)
anova(m1, m2)

# Likelihood ratio tests of Negative Binomial Models
# 
# Response: ent_count
# Model     theta Resid. df    2 x log-lik.   Test    df
# 1      grass.nograss 0.5568987       674       -6227.537             
# 2 grass.nograss + do 0.5878148       673       -6184.997 1 vs 2     1
# LR stat.      Pr(Chi)
# 1                      
# 2 42.53946 6.927059e-11

###### Checking model assumptions ####
# estimating a dispersion parameter (not shown in the output) 
# that is held constant in a Poisson model
m3 <- glm(ent_count ~ grass.nograss + do, family = "poisson", data = dat)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)

## 'log Lik.' 0 (df=4)


## get the confidence intervals for the coefficients 
## by profiling the likelihood function
(est <- cbind(Estimate = coef(m1), confint(m1)))
# Waiting for profiling to be done...
# Estimate       2.5 %     97.5 %
#   (Intercept)     5.3781619  4.85467926  5.8951775
# grass.nograss1  0.2638860  0.05965523  0.4681696
# do             -0.2319473 -0.29654381 -0.1655627

# Look at incident rate ratios 
# exponentiate our model coefficients. 
# same applies to the confidence intervals.

exp(est)

# Estimate       2.5 %      97.5 %
#   (Intercept)    216.6237434 128.3395210 363.2813249
# grass.nograss1   1.3019798   1.0614705   1.5970683
# do               0.7929879   0.7433831   0.8474167

# predicted values - can't do this with site and grass.nograss

#run
newdata1 <- data.frame(do = mean(dat$do), 
                       grass.nograss = factor(1:2, levels = 1:2, 
                                     labels = levels(dat$grass.nograss)))
newdata1$phat <- predict(m1, newdata1, type = "response")
newdata1

# do grass.nograss     phat
# 1 8.030991             0 33.62922
# 2 8.030991             1 43.78456

## obtain the mean predicted number of events for values of grass no grass
# across its entire range for each level of site and graph these

newdata2 <- data.frame(
  do = rep(seq(from = min(dat$do), to = max(dat$do), 
                 length.out = 96), 3),
  grass.nograss = factor(rep(1:2, each = 100), levels = 1:2, labels =
                  levels(dat$grass.nograss)))

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