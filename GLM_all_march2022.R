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

# plot the data
ggplot(dat, aes(ent_count, fill = dist.from.shore)) + geom_histogram(binwidth = 1) +
  facet_grid(dist.from.shore ~ ., margins = TRUE, scales = "free")

#explore data
with(dat, tapply(grass.nograss, grass, nograss))

#GLM analysis 
summary(m1 <- glm.nb(ent_count ~ grass.nograss + site + season + dist.from.shore + do,
                     data = dat))
# Call:
#   glm.nb(formula = ent_count ~ grass.nograss + site + season + 
#            dist.from.shore + do, data = dat, init.theta = 1.121919487, 
#          link = log)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -3.2950  -1.0264  -0.3363   0.3874   4.0471  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         4.25936    0.28148  15.132  < 2e-16 ***
#   grass.nograss1      0.17244    0.07582   2.274 0.022946 *  
#   siteGG              2.06776    0.10458  19.771  < 2e-16 ***
#   sitesmithcove       2.16107    0.10597  20.393  < 2e-16 ***
#   seasonspring       -0.42385    0.17261  -2.456 0.014066 *  
#   seasonsummer       -1.03601    0.11641  -8.900  < 2e-16 ***
#   seasonwinter       -0.89799    0.11796  -7.613 2.68e-14 ***
#   dist.from.shore20  -0.31591    0.12624  -2.503 0.012331 *  
#   dist.from.shore40  -0.60340    0.12687  -4.756 1.97e-06 ***
#   dist.from.shore60  -0.48229    0.12717  -3.793 0.000149 ***
#   dist.from.shore80  -0.61121    0.12715  -4.807 1.53e-06 ***
#   dist.from.shore100 -0.65607    0.12880  -5.094 3.51e-07 ***
#   do                 -0.19722    0.03668  -5.377 7.58e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(1.1219) family taken to be 1)
# 
# Null deviance: 1527.8  on 711  degrees of freedom
# Residual deviance:  825.2  on 699  degrees of freedom
# (116 observations deleted due to missingness)
# AIC: 6095.9
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  1.1219 
# Std. Err.:  0.0627 
# 
# 2 x log-likelihood:  -6067.9020   

# see if site is significant
m2 <- update(m1, . ~ . - site)
anova(m1, m2)

# Likelihood ratio tests of Negative Binomial Models
# 
# Response: ent_count
# Model     theta Resid. df
# 1        grass.nograss + season + dist.from.shore + do 0.7469141       701
# 2 grass.nograss + site + season + dist.from.shore + do 1.1219195       699
# 2 x log-lik.   Test    df LR stat. Pr(Chi)
# 1       -6371.402                              
# 2       -6067.902 1 vs 2     2 303.4995       0 (site increases log likelihood, good)

## season
m2 <- update(m1, . ~ . - season)
anova(m1, m2)

# Likelihood ratio tests of Negative Binomial Models
# 
# Response: ent_count
# Model     theta Resid. df
# 1          grass.nograss + site + dist.from.shore + do 0.9689313       702
# 2 grass.nograss + site + season + dist.from.shore + do 1.1219195       699
# 2 x log-lik.   Test    df LR stat. Pr(Chi)
# 1       -6168.344                              
# 2       -6067.902 1 vs 2     3  100.442       0 (season increases log likelihood, good)

## dist from shore
m2 <- update(m1, . ~ . - dist.from.shore)
anova(m1, m2)
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: ent_count
# Model    theta Resid. df    2 x log-lik.
# 1                   grass.nograss + site + season + do 1.064784       704       -6106.299
# 2 grass.nograss + site + season + dist.from.shore + do 1.121919       699       -6067.902
# Test    df LR stat.      Pr(Chi)
# 1                                   
# 2 1 vs 2     5 38.39664 3.140549e-07 (dist from shore increases log likelihood, good)

## grass no grass
m2 <- update(m1, . ~ . - grass.nograss)
anova(m1, m2)

# Likelihood ratio tests of Negative Binomial Models
# 
# Response: ent_count
# Model    theta Resid. df    2 x log-lik.
# 1                 site + season + dist.from.shore + do 1.114059       700       -6073.083
# 2 grass.nograss + site + season + dist.from.shore + do 1.121919       699       -6067.902
# Test    df LR stat.    Pr(Chi)
# 1                                 
# 2 1 vs 2     1 5.180528 0.02284138 (grass.nograss increases log likelihood, good)

## DO
m2 <- update(m1, . ~ . - do)
anova(m1, m2)
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: ent_count
# Model     theta Resid. df
# 1      grass.nograss + site + season + dist.from.shore 0.9620343       795
# 2 grass.nograss + site + season + dist.from.shore + do 1.1219195       699
# 2 x log-lik.   Test    df LR stat. Pr(Chi)
# 1       -6669.229                              
# 2       -6067.902 1 vs 2    96  601.327       0 (Do increases log likelihood, good)

###### Checking model assumptions ####
# estimating a dispersion parameter (not shown in the output) 
# that is held constant in a Poisson model
m3 <- glm(ent_count ~ grass.nograss + site + season + dist.from.shore, 
          family = "poisson", data = dat)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)

## 'log Lik.' 0 (df=14)


## get the confidence intervals for the coefficients 
## by profiling the likelihood function
(est <- cbind(Estimate = coef(m1), confint(m1)))
# Waiting for profiling to be done...
#                      Estimate       2.5 %      97.5 %
#   (Intercept)         4.2593637  3.68734362  4.82986958
# grass.nograss1      0.1724435  0.02419158  0.32072741
# siteGG              2.0677572  1.86203994  2.27179153
# sitesmithcove       2.1610739  1.93902644  2.38298359
# seasonspring       -0.4238548 -0.78767982 -0.05895258
# seasonsummer       -1.0360117 -1.28111857 -0.79109962
# seasonwinter       -0.8979898 -1.12253533 -0.67564678
# dist.from.shore20  -0.3159109 -0.56541022 -0.06658260
# dist.from.shore40  -0.6034026 -0.85422236 -0.35280224
# dist.from.shore60  -0.4822945 -0.73479448 -0.22990268
# dist.from.shore80  -0.6112105 -0.86463569 -0.35794144
# dist.from.shore100 -0.6560661 -0.91311770 -0.39891043
# do                 -0.1972199 -0.27141940 -0.12202543

# Look at incident rate ratios 
# exponentiate our model coefficients. 
# same applies to the confidence intervals.

exp(est)

#                      Estimate      2.5 %      97.5 %
#   (Intercept)        70.7649415 39.9386138 125.1946315
# grass.nograss1      1.1882046  1.0244866   1.3781299
# siteGG              7.9070688  6.4368542   9.6967572
# sitesmithcove       8.6804550  6.9519795  10.8371884
# seasonspring        0.6545189  0.4548990   0.9427515
# seasonsummer        0.3548672  0.2777265   0.4533460
# seasonwinter        0.4073878  0.3254536   0.5088272
# dist.from.shore20   0.7291244  0.5681270   0.9355856
# dist.from.shore40   0.5469474  0.4256140   0.7027162
# dist.from.shore60   0.6173652  0.4796040   0.7946109
# dist.from.shore80   0.5426936  0.4212050   0.6991140
# dist.from.shore100  0.5188886  0.4012712   0.6710508
# do                  0.8210101  0.7622967   0.8851259

# predicted values - can't do this with site and grass.nograss
# make site a factor
dat$site <- as.factor(dat$site)
# make dist from shore a factor
dat$dist.from.shore <- as.factor(dat$dist.from.shore)
# make grass a factor
dat$grass.nograss <- as.factor(dat$grass.nograss)
# make season from shore a factor
dat$season <- as.factor(dat$season)
#run
newdata1 <- data.frame(dist.from.shore = factor(1:6, levels = 1:6, 
                                                labels = levels(dat$dist.from.shore)),
                       grass.nograss = factor(1:2, levels = 1:2, 
                                              labels = levels(dat$grass.nograss)),
                       season = factor(1:4, levels = 1:4, 
                                       labels = levels(dat$season)),
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