# 2/24/20
# field data

#clear workspace
rm(list=ls())

# call in data 
# from downloads
data <- read.csv("Field Assay Counts - Sheet1.csv")
# view data
head(data)
#remove geom data 
data <- data[-c(6,7,8,9,10,14,15)]
#plot temp data
plot(data$Date, data$temp)
#remove discoveryfield site data
data <- data[data$site != "Disc",]
# create data for each site
Alki <- data[data$site=='Alki',]
Alki_ent <- na.omit(Alki$ent_count)
GG <- data[data$site=='GG',]
GG_ent <- na.omit(GG$ent_count)
Smith <- data[data$site=="Smith's Cove",]
Smith_ent <- na.omit(Smith$ent_count)
# bacteria per site
# averages per site
a_avg <- mean(Alki_ent)
g_avg <- mean(GG_ent)
s_avg <- mean(Smith_ent)


# bar plot
library(dplyr)
library(ggplot2)
#create sd 
#sds <- c(sd(Alki_ent), sd(GG_ent), sd(Smith_ent))
#means <- c(a_avg, g_avg, s_avg)
#name <- c("Alki","GG","Smith's Cove")
#errors <- data.frame(site = name, mean = means, lower = sds - means, 
                    #upper = sds + means)
#remove NAs
data <- na.omit(data)
#summarize data
data_summary <- data %>% # the names of the new data frame and the data frame to be summarised
   group_by(site) %>% # the grouping variable
   summarise(
      mean_ent = mean(ent_count), # calculates the mean of each group
      sd_ent = sd(ent_count), # calculates the standard deviation of each group
      n_ent = n(), # calculates the sample size per group
      SE_ent = sd(ent_count) / sqrt(n())
   ) # calculates the standard error of each group
#plot data
EntPlot <- ggplot(data_summary, aes(site, mean_ent)) +
   geom_col() +
   geom_errorbar(aes(ymin = mean_ent - sd_ent, ymax = mean_ent + sd_ent), width = 0.2)

EntPlot + labs(y = "Enterococcus abundance (CFU/100ml)", x = "Site") + theme_classic()


#create new data frame with just grass and counts
grass <- data[which(data$ent_count & data$grass.nograss == 1), ]
ent <- grass$ent_count
mean(ent) #41.05
sd(ent)#57.01
length(grass$ent_count) #390

no_grass <- data[which(data$ent_count & data$grass.nograss == 0), ]
ent2 <- no_grass$ent_count
mean(ent2) #37.76
sd(ent2) #44.80
length(no_grass$ent_count) #359
# plot ent count by dist from shore

#PLOT
# grass present
plot(grass$sample, grass$ent_count, main = "Ent counts by Dist from Shore,eelgrass present",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,500))
abline(h=104, col ="red")
#grass absent 
plot(no_grass$sample, no_grass$ent_count, main = "Ent counts by Dist from Shore,no eelgrass present",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,500))
abline(h=104, col ="red")

#stats 2
t.test(x=grass$ent_count,y=no_grass$ent_count, alternative = "less",
       mu=0,paired = FALSE, var.equal = T, conf.level = 0.95)

#ttest
t =(22.099-22.835)/(29.666/sqrt(202)) #-0.0017456

#plot data 
# ent counts by transect 
plot(data$transect, data$ent_count, 
     xlab = "Transect", ylab = "Enterococcus Single Sample Count",
     main = "Ent. Counts by transect")
# ent counts by treatment
plot(data$grass.nograss, data$ent_count, 
     xlab = "Seagrass coverage", ylab = "Enterococcus Single Sample Count (Cfu/100ml)",
     main = "Ent. Counts by treatment")
# ent counts by site
plot(data$site, data$ent_count, 
     xlab = "Site Name", ylab = "Enterococcus Single Sample Count (Cfu/100ml)",
     main = "Enterococcus Counts by Site")
abline(h=104, col ="red")
#ent count by season
plot(data$season, data$ent_count, 
     xlab = "Season", ylab = "Enterococcus Single Sample Count (Cfu/100ml)",
     main = "Enterococcus Counts by Season")
abline(h=104, col ="red")
# geometric mean by site 
plot(data$transect, data$geom_mean, 
     xlab = "Transect", ylab = "Enterococcus geom mean",
     main = "Ent. Counts by transect")
# geometric mean by site
plot(data$site, data$geom_mean, 
     xlab = "Site", ylab = "Enterococcus geom mean",
     main = "Ent. Counts by Site")


#geom mean by grass/no grass
plot(data$grass.nograss, data$geom_mean, 
     xlab = "Seagrass coverage", ylab = "Enterococcus geom mean",
     main = "Ent. Counts by treatment")
#remove NA's
newdat <- data[217-252,]
plot(newdat$Date, newdat$ent_count,
     xlab = "Collection Date", ylab = "Enterococcus Single Sample Count (Cfu/100ml)",
     main = "Enterococcus Counts by Date")
abline(h=104, col ="red")
#plot by season
plot(newdat$season, newdat$ent_count,
     xlab = "Collection Date", ylab = "Enterococcus Single Sample Count (Cfu/100ml)",
     main = "Enterococcus Counts by Date")
abline(h=104, col ="red")
#plot by site
nodisc <- newdat[which(newdat$site != "Disc"),]
#plot ent by site 
nodisc$site <- as.factor(nodisc$site) # need to code site as a factor!!!
plot(nodisc$site, nodisc$ent_count,
     xlab = "Site", ylab = "Enterococcus Single Sample Count (Cfu/100ml)",
     main = "Enterococcus Counts by Site")
abline(h=104, col ="red")
# salinity proxy of TDS?
plot(nodisc$grass.nograss, nodisc$sal)
# data per site 
Alki <- newdat[which(newdat$site=='Alki'),]
GG <- newdat[which(newdat$site=='GG'),]
Disc <- newdat[which(newdat$site=='Disc'),]
Smith <- newdat[which(newdat$site=="Smith's Cove"),]

#plot temp data for each site
#corellations 
#alki 
plot(Alki$temp, Alki$ent_count, main ="Alki Temp vs Ent Count",
     xlab = "Temperature", 
     ylab = "Ent Count")
plot(Alki$ph, Alki$ent_count, main ="Alki pH vs Ent Count",
     xlab = "pH", 
     ylab = "Ent Count")
plot(Alki$do, Alki$ent_count, main ="Alki DO vs Ent Count",
     xlab = "Dissolved Oxygen", 
     ylab = "Ent Count")
plot(Alki$par, Alki$ent_count, main ="Alki PAR vs Ent Count",
     xlab = "PAR", 
     ylab = "Ent Count")
plot(Alki$sal, Alki$ent_count, main ="Alki Salinity vs Ent Count",
     xlab = "Salinity", 
     ylab = "Ent Count")
#GG 
plot(GG$temp, GG$ent_count, main ="GG Temp vs Ent Count",
     xlab = "Temperature", 
     ylab = "Ent Count")
plot(GG$ph, GG$ent_count, main ="GG pH vs Ent Count",
     xlab = "pH", 
     ylab = "Ent Count")
plot(GG$do, GG$ent_count, main ="GG DO vs Ent Count",
     xlab = "Dissolved Oxygen", 
     ylab = "Ent Count")
plot(GG$par, GG$ent_count, main ="GG PAR vs Ent Count",
     xlab = "PAR", 
     ylab = "Ent Count")
plot(GG$sal, GG$ent_count, main ="GG Salinity vs Ent Count",
     xlab = "Salinity", 
     ylab = "Ent Count")
#sc
plot(Smith$temp, Smith$ent_count, main ="SC Temp vs Ent Count",
     xlab = "Temperature", 
     ylab = "Ent Count")
plot(Smith$ph, Smith$ent_count, main ="SC pH vs Ent Count",
     xlab = "pH", 
     ylab = "Ent Count")
plot(Smith$do, Smith$ent_count, main ="SC DO vs Ent Count",
     xlab = "Dissolved Oxygen", 
     ylab = "Ent Count")
plot(Smith$par, Smith$ent_count, main ="SC PAR vs Ent Count",
     xlab = "PAR", 
     ylab = "Ent Count")
plot(Smith$sal, Smith$ent_count, main ="SC Salinity vs Ent Count",
     xlab = "Salinity", 
     ylab = "Ent Count")
#water chem thru time 
#temp
plot(Alki$Date, Alki$temp, main = "Alki Temp Through Time",
          xlab = "Date", 
          ylab = "Temperature")
plot(GG$Date, GG$temp, main = "GG Temp Through Time",
     xlab = "Date", 
     ylab = "Temperature")
plot(Smith$Date, Smith$temp, main = "Smith's Cove Temp Through Time",
     xlab = "Date", 
     ylab = "Temperature")
# pH
plot(Alki$Date, Alki$pH, main = "Alki Temp Through Time",
     xlab = "Date", 
     ylab = "pH")
plot(GG$Date, GG$pH, main = "GG Temp Through Time",
     xlab = "Date", 
     ylab = "pH")
plot(Smith$Date, Smith$pH, main = "Smith's Cove Temp Through Time",
     xlab = "Date", 
     ylab = "pH")
# DO
plot(Alki$Date, Alki$do, main = "Alki Temp Through Time",
     xlab = "Date", 
     ylab = "pH")
plot(GG$Date, GG$do, main = "GG Temp Through Time",
     xlab = "Date", 
     ylab = "pH")
plot(Smith$Date, Smith$do, main = "Smith's Cove Temp Through Time",
     xlab = "Date", 
     ylab = "pH")

# site change over time
plot(Alki$Date, Alki$ent_count, type="l", col="blue",main = "Alki Beach Park")
plot(GG$Date, GG$ent_count, type="l", col="red",main = "Golden Gardens Park")
plot(Disc$Date, Disc$ent_count, type="l", col="green",main = "Discovery Park")
plot(Smith$Date, Smith$ent_count, type="l", col="grey",main = "Smith's Cove",
     xlab = "Date", ylab = "Enterococcus Single Sample Count (Cfu/100ml)")

#count by dist from shore
plot(Alki$sample, Alki$ent_count, main = "Alki Beach Park",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
plot(GG$sample, GG$ent_count, main = "Alki Beach Park",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
#site counts by transect
plot(Alki$transect, Alki$ent_count, main = "Alki Beach Park",
     xlab = "Transect", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,150))
plot(GG$transect, GG$ent_count, main = "Golden Gardens Park",
     xlab = "Transect", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)",ylim = c(0,150))
plot(Disc$transect, Disc$ent_count, main = "Discovery Park",
     xlab = "Transect", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,150))
plot(Smith$transect, Smith$ent_count, main = "Smith's Cove",
     xlab = "Transect", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)",ylim = c(0,150))

#plot grass vs. no grass
plot(Alki$grass.nograss, Alki$ent_count, main = "Alki Beach Park",
     xlab = "Treatment", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")

plot(GG$grass.nograss, GG$ent_count, main = "GG",
     xlab = "Treatment", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
plot(Disc$grass.nograss, Disc$ent_count, main = "Discovery",
     xlab = "Treatment", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
plot(Smith$grass.nograss, Smith$ent_count, main = "Smith's Cove",
     xlab = "Treatment", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,200))

abline(h=104, col ="red")
plot(newdat$grass.nograss, newdat$ent_count, main = "Overall",
      xlab = "Treatment", 
      ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
#counts over time by site
#Alki
plot(Alki$Date[which(Alki$grass.nograss == 1)], # grass present
     Alki$ent_count[which(Alki$grass.nograss == 1)], 
     main = "Alki Beach Park, Eelgrass Present",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
plot(Alki$Date[which(Alki$grass.nograss == 0)], # grass absent
     Alki$ent_count[which(Alki$grass.nograss == 0)], 
     main = "Alki Beach Park, Eelgrass Absent",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")

plot(Alki$Date, Alki$ent_count, 
     main = "Alki Beach Park",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,300))
abline(h=104, col ="red")
#gg
plot(GG$Date[which(GG$grass.nograss == 1)], # grass present
     GG$ent_count[which(GG$grass.nograss == 1)], 
     main = "GG Park, Eelgrass Present",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,200))
abline(h=104, col ="red")
plot(GG$Date[which(GG$grass.nograss == 0)], # grass absent
     GG$ent_count[which(GG$grass.nograss == 0)], 
     main = "GG Park, Eelgrass Absent",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,200))
abline(h=104, col ="red")

plot(GG$Date, 
     GG$ent_count, 
     main = "Golden Garden Park",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,300))
abline(h=104, col ="red")
#Smith's Cove
plot(Smith$Date[which(Smith$grass.nograss == 1)], # grass present
     Smith$ent_count[which(Smith$grass.nograss == 1)], 
     main = "Smith's Cove, Eelgrass Present",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,300))
abline(h=104, col ="red")
plot(Smith$Date[which(Smith$grass.nograss == 0)], # grass absent
     Smith$ent_count[which(Smith$grass.nograss == 0)], 
     main = "Smith's Cove, Eelgrass Absent",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,300))
abline(h=104, col ="red")

plot(Smith$Date, 
     Smith$ent_count, 
     main = "Smith's Cove",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,300))
abline(h=104, col ="red")
# Discovery
plot(Disc$Date[which(Disc$grass.nograss == 1)], # grass present
     Disc$ent_count[which(Disc$grass.nograss == 1)], 
     main = "Discovery Park, Eelgrass Present",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
plot(Disc$Date[which(Disc$grass.nograss == 0)], # grass absent
     Disc$ent_count[which(Disc$grass.nograss == 0)], 
     main = "Discovery Park, Eelgrass Absent",
     xlab = "Date", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
#check water chem 
plot(newdat$grass.nograss, newdat$temp, main = "temp comparison")
plot(newdat$grass.nograss, newdat$ph, main = "pH comparison")
plot(newdat$grass.nograss, newdat$sal, main = "Salinity comparison")
plot(newdat$grass.nograss, newdat$do, main = "D.O. comparison")

plot(newdat$Date, newdat$temp)

#site stats
#Alki
alki2 <- Alki[which(Alki$ent_count & Alki$grass.nograss == 1), ]
A_ent <- alki2$ent_count
mean(A_ent) # 5.17
sd(A_ent)#5.42
length(alki2$ent_count) #46
alki3 <- Alki[which(Alki$ent_count & Alki$grass.nograss == 0), ]
B_ent <- alki3$ent_count
mean(B_ent) # 5.83



#count by dist from shore
plot(alki2$sample, alki2$ent_count, main = "Alki Beach Park eelgrass",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")

plot(alki3$sample, alki3$ent_count, main = "Alki Beach Park no eelgrass",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")


#histogram
hist(alki2$ent_count)
hist(alki3$ent_count)

#ttest
t =(5.17-5.83)/(5.42/sqrt(46)) #-.01795918
pval <- pt(q = -.82589, df = 45) #.2066, not sig
#t.test(x=alki2$ent_count,y=alki3$ent_count, alternative = "less",
       #mu=0,paired = FALSE, var.equal = T, conf.level = 0.95)
#disc
disc2 <- Disc[which(Disc$ent_count & Disc$grass.nograss == 1), ]
D_ent <- disc2$ent_count
mean(D_ent) #5.75
sd(D_ent)#7.02
length(disc2$ent_count) #16
disc3 <- Disc[which(Disc$ent_count & Disc$grass.nograss == 0), ]
E_ent <- disc3$ent_count
mean(E_ent) #11.25
#count by dist from shore
plot(disc2$sample, disc2$ent_count, main = "Discovery Park eelgrass",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")

plot(disc3$sample, disc3$ent_count, main = "Discovery Park no eelgrass",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
#histogram
hist(disc2$ent_count)
hist(disc3$ent_count)
#ttest
t =(5.75-11.25)/(7.02/sqrt(16)) #-.1958689
pt(q = -3.1339, df = 15) #.0034 (sig!!)
#t.test(x=disc2$ent_count,y=disc3$ent_count, alternative = "less",
       #mu=0,paired = FALSE, var.equal = T, conf.level = 0.95)


#GG
GG2 <- GG[which(GG$ent_count & GG$grass.nograss == 1), ]
G_ent <- GG2$ent_count
mean(G_ent) #42.59722
sd(G_ent)#40.42091
length(GG2$ent_count) #72
GG3 <- GG[which(GG$ent_count & GG$grass.nograss == 0), ]
H_ent <- GG3$ent_count
mean(H_ent) #37.08826
#count by dist from shore
plot(GG2$sample, GG2$ent_count, main = "GG Park eelgrass",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")

plot(GG3$sample, GG3$ent_count, main = "GG Park no eelgrass",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
#histogram
hist(GG2$ent_count)
hist(GG3$ent_count)
#ttest
t =(42.59722-37.08826)/(40.42091/sqrt(72)) #.01606191
pt(q = 1.1565, df = 71) #.87 (not sig)
t.test(x=GG2$ent_count,y=GG3$ent_count, alternative = "less",
       mu=0,paired = FALSE, var.equal = T, conf.level = 0.95)
#Smith's Cove
SC2 <- Smith[which(Smith$ent_count & Smith$grass.nograss == 1), ]
SC_ent <- SC2$ent_count
mean(SC_ent) #15.69118
sd(SC_ent)#11.15061
length(SC2$ent_count) #68
SC3 <- Smith[which(Smith$ent_count & Smith$grass.nograss == 0), ]
S_ent <- SC3$ent_count
mean(S_ent) #19.54717
#count by dist from shore
plot(SC2$sample, SC2$ent_count, main = "GG Park eelgrass",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")

plot(SC3$sample, SC3$ent_count, main = "GG Park no eelgrass",
     xlab = "Dist from Shore", 
     ylab = "Enterococcus Single Sample Count (Cfu/100ml)", ylim = c(0,160))
abline(h=104, col ="red")
#histogram
hist(SC2$ent_count)
hist(SC3$ent_count)
#ttest
t =(15.69118-19.54717)/(11.15061/sqrt(68)) #-.0419356
pt(q = -2.85, df = 67) #.0029 (sig!!)
#t.test(x=SC2$ent_count,y=SC3$ent_count, alternative = "less",
       #mu=19.5,paired = FALSE, var.equal = T, conf.level = 0.95)


