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
ggplot(dat, aes(ent_count, fill = season)) + geom_histogram(binwidth = 1) +
  facet_grid(season ~ ., margins = TRUE, scales = "free")

#GLM analysis 
summary(m1 <- glm.nb(ent_count ~ site + season, data = dat))
## Call:
# glm.nb(formula = ent_count ~ par + do, data = dat, init.theta = 0.6095931712, 
#        link = log)

# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3594  -1.1635  -0.4039   0.2236   3.0564  

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  5.021659   0.282403  17.782  < 2e-16 ***
#   par         -0.010652   0.002009  -5.301 1.15e-07 ***
#   do          -0.150631   0.035816  -4.206 2.60e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for Negative Binomial(0.6096) family taken to be 1)

# Null deviance: 883.58  on 675  degrees of freedom
# Residual deviance: 808.26  on 673  degrees of freedom
# (152 observations deleted due to missingness)
# # AIC: 6165.1

# Number of Fisher Scoring iterations: 1


# Theta:  0.6096 
# Std. Err.:  0.0313 

# 2 x log-likelihood:  -6157.0920 Call:


# see if PAR is significant
m2 <- update(m1, . ~ . - par)
anova(m1, m2)

# Likelihood ratio tests of Negative Binomial Models

# Response: ent_count
# Model     theta Resid. df    2 x log-lik.   Test    df LR stat. Pr(Chi)
# 1       do 0.6125512       710       -6527.827                              
# 2 par + do 0.6095932       673       -6157.092 1 vs 2    37 370.7352       0


## see is DO is signficant 
m2 <- update(m1, . ~ . - do)
anova(m1, m2)
## Likelihood ratio tests of Negative Binomial Models

## Response: ent_count
#  Model     theta Resid. df    2 x log-lik.   Test    df LR stat. Pr(Chi)
#  1      par 0.5272816       733       -6541.034                              
#  2 par + do 0.6095932       673       -6157.092 1 vs 2    60 383.9423       0

###### Checking model assumptions ####
# estimating a dispersion parameter (not shown in the output) 
# that is held constant in a Poisson model
m3 <- glm(ent_count ~ par + do, family = "poisson", data = dat)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)

## 'log Lik.' 0 (df=4)


## get the confidence intervals for the coefficients 
## by profiling the likelihood function
(est <- cbind(Estimate = coef(m1), confint(m1)))
## Waiting for profiling to be done...
## Estimate       2.5 %       97.5 %
## (Intercept)  5.02165863  4.49955900  5.536542413
## par         -0.01065167 -0.01356841 -0.007466146
## do          -0.15063097 -0.21482394 -0.084456259

# Look at incident rate ratios 
# exponentiate our model coefficients. 
# same applies to the confidence intervals.

exp(est)

##             Estimate      2.5 %      97.5 %
## (Intercept) 151.6626478 89.9774427 253.7989486
## par           0.9894049  0.9865232   0.9925617
## do            0.8601651  0.8066835   0.9190119

# predicted values - 
# make site a factor
dat$site <- as.factor(dat$site)

#run
newdata1 <- data.frame(par = mean(dat$par), 
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