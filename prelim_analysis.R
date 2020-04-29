##These are the SI analyses of diet of mink, marten, weasel, and otter in JPRF ----
#performed by D. Breault in 2019-2020 for Dexter - last updated 14apr2020

#set library

.libPaths("C:/Google Drive/R/library")

#set working directory

setwd(("C:/Google Drive/R/projects/jprf_si_analysis"))

#remove any objects from previous analysis

rm(list = ls())

#Install packages ----

install.packages("devtools")

install.packages("extrafont")

install.packages("ggpubr")

install.packages("PMCMR")

install.packages("sqldf")

install.packages('vegan')

install.packages("siar")

install.packages("nicheROVER")

install.packages("dplyr")

install.packages('ggstance')

install.packages('ggvis')

install.packages('class')

install.packages('gmodels')

install.packages("data.table")

install.packages("R")

install.packages("ggdistribute")

install.packages("ggproto")

install.packages("hrbrthemes")

install.packages("viridis")

install.packages("sjPlot")

install.packages("stargazer")

#prep (load packages) ----

#Load packages


library("ggplot2")

library("hrbrthemes")

library("tidyr")

library("viridis")

library('devtools')

library('tidyverse')

library('lubridate')

library('ggpubr')

library('siar')

library('nicheROVER')

library('plyr')

library('dplyr')

library('reshape2')

library('PMCMR')

library('sqldf')

library('vegan')

library('ggstance')

library('class')

library('gmodels')

library('data.table')

library('R')

library("ggdistribute")

library("sjPlot")

library("stargazer")

library("extrafont")

font_import()

y
n

theme_set(theme_pubr())

#Register fonts for Windows bitmap output

loadfonts(device="win")      

#Is hair significantly different in C & N btwn species?----
#I did this on December 8th, 2019-April, 28th, 2020

#load SI data (cleaned) from hair of consumers
#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                     sql = "select * from file where `Tissue` = 'hair'", sep = ",")

summary(dat.hair)

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

boxplot(dat.hair$C)

boxplot(dat.hair$N)

# Shapiro-Wilk normality test for C13 and N15 by species

with(dat.hair, shapiro.test(N[spps_code == 1]))

with(dat.hair, shapiro.test(C[spps_code == 1]))

with(dat.hair, shapiro.test(N[spps_code == 2]))

with(dat.hair, shapiro.test(C[spps_code == 2]))

with(dat.hair, shapiro.test(N[spps_code == 3]))

with(dat.hair, shapiro.test(C[spps_code == 3]))

with(dat.hair, shapiro.test(N[spps_code == 4]))

with(dat.hair, shapiro.test(C[spps_code == 4]))

### marten SI data do not follow a normal distribution (p < 0.01)

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(C ~ spps_code, data = dat.hair)

require(PMCMR)

data(dat.hair)

attach(dat.hair)

posthoc.kruskal.nemenyi.test(x = C, g = spps_code, dist="Tukey")

###no diff in C btwn mink and marten (p = .1), but diff btwn mink & weasel (p << 0.001)
###       and btwn marten and weasel (p = .02) and btwn weasel and otter (p = 0.004)
##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(N ~ spps_code, data = dat.hair)

require(PMCMR)

data(dat.hair)

attach(dat.hair)

posthoc.kruskal.nemenyi.test(x = N, g = spps_code, dist="Tukey")

#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

###diff btwn mink and marten (p << 0.001) btwn mink and weasel (p = 0.01)
###btwn marten and otter (p << 0.001) and btwn weasel and otter (p = 0.002)

#remove any objects from previous analysis

rm(list = ls())

#plot the mustelid hair data----

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

summary(dat.hair)

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

boxplot(dat.hair$C)

boxplot(dat.hair$N)

# Join means by species to the original data frame and pipe to ggplot

consumer_means <- left_join(dat.hair, 
                    dat.hair %>%
                     group_by(species) %>%
                     summarise_at(vars(C, N), funs(mean = mean))
                ) %>% 
                      ggplot(aes(colour=species, shape = species)) +
                      # Plot raw data
                      geom_point(aes(x=C, y=N), size = 2, alpha = 0.7) +
                      # Plot mean by species
                      geom_point(aes(x=C_mean, y=N_mean), size=3) +
                      # Add vertical errorbars 
                      scale_shape_manual(values = c(0, 21, 17, 18)) +
                      scale_color_manual(values = c("gray1", "gray1", "gray15", "gray1")) + 
                      stat_summary(fun.data = mean_se, fun.args = list(mult=1.96), 
                                   aes(x = C_mean, y=N),
                                   geom="errorbar", width=0.05) +
  # Add horizontal errorbars
  stat_summaryh(fun.data=mean_se_h, fun.args=list(mult=1.96), 
                aes(x=C, y=N_mean),
                geom="errorbarh", width=0.1) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Isotopic Signatures", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "16",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.88, 0.98)) +
  scale_x_continuous(breaks = seq(-30, -20, 1)) +
  scale_y_continuous(breaks = seq(0, 14, 1))
consumer_means <- consumer_means + labs(colour = "Mustelid Species", shape = "Mustelid Species")

  tiff(file = "figures/consumer_mean_95CI_18apr2020.tif", 
     units = "cm", width = 20, height = 15, res = 300)

consumer_means

dev.off()

#Niche Ellipses with mustelid species as consumer groupings----

# load SI data (cleaned) from hair of consumers
#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

outlier <- boxplot(dat.hair$N, plot = FALSE)$out

dat.hair[which(dat.hair$N %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$N %in% outlier),]

boxplot(dat.hair$C)

boxplot(dat.hair$N)

#summarize by species

summary.mustelid <- as.data.frame(dat.hair %>%
                                group_by(species) %>%
                                dplyr::summarize(mean_C = mean(C),
                                                 se_C = sd(C)/sqrt(n()),
                                                 mean_N = mean(N),
                                                 se_N=sd(N)/sqrt(n()),
                                                 sample_size = n()))

summary.mustelid

write.csv(summary.mustelid, "outputs/summary_mustelids_28apr2020.csv")

#Species Codes
#mink = 1, marten = 2, weasel = 3, otter = 4

#use mustelid data to calculate standard ellipse metrics (SEAc) by species

mink <- subset(dat.hair, species == 'mink')

mink_SEA <- data.frame(standard.ellipse(mink$C, mink$N, confs = NULL, steps = 5))

write.csv(mink_SEA, "outputs/standard.ellipse statistics/mink_SEA_28apr2020.csv")

f_mink <- subset(dat.hair, species == 'mink' & sex == 'F')

f_mink_SEA <- data.frame(standard.ellipse(f_mink$C, f_mink$N, confs = NULL, steps = 5))

write.csv(f_mink_SEA, "outputs/standard.ellipse statistics/f_mink_SEA_28apr2020.csv")

m_mink <- subset(dat.hair, species == 'mink' & sex == 'M')

m_mink_SEA <- data.frame(standard.ellipse(m_mink$C, m_mink$N, confs = NULL, steps = 5))

write.csv(m_mink_SEA, "outputs/standard.ellipse statistics/m_mink_SEA_28apr2020.csv")

marten <- subset(dat.hair, spps_code == '2')

marten_SEA <- data.frame(standard.ellipse(marten$C, marten$N, confs = NULL, steps = 5))

write.csv(marten_SEA, "outputs/standard.ellipse statistics/marten_SEA_28apr2020.csv")

weasel <- subset(dat.hair, spps_code == '3')

weasel_SEA <- data.frame(standard.ellipse(weasel$C, weasel$N, confs = NULL, steps = 5))

write.csv(weasel_SEA, "outputs/standard.ellipse statistics/weasel_SEA_28apr2020.csv")

otter <- subset(dat.hair, spps_code == '4')

otter_SEA <- data.frame(standard.ellipse(otter$C, otter$N, confs = NULL, steps = 5))

write.csv(otter_SEA, "outputs/standard.ellipse statistics/otter_SEA_28apr2020.csv")

#remove any objects from previous analysis

rm(list = ls())

#plot the SEAc for each species----

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#Find and remove outliers from data set

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

outlier <- boxplot(dat.hair$N, plot = FALSE)$out

dat.hair[which(dat.hair$N %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$N %in% outlier),]

#subset consumer data by species manually

mink <- subset(dat.hair, species == 'mink')

marten <- subset(dat.hair, species == 'marten')

weasel <- subset(dat.hair, species == 'weasel')

otter <- subset(dat.hair, species == 'otter')

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(mink$C, mink$N)
SE2 <- standard.ellipse(marten$C, marten$N)
SE3 <- standard.ellipse(weasel$C, weasel$N)
SE4 <- standard.ellipse(otter$C, otter$N)

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

mink_ <- rep("mink", length(SE1$xSEAc))
marten_ <- rep("marten", length(SE2$xSEAc))
weasel_ <- rep("weasel", length(SE3$xSEAc))
otter_ <- rep("otter", length(SE4$xSEAc))

#create new data frame with names and ellipse outputs

species <- c(mink_,marten_,weasel_,otter_)
x <- c(SE1$xSEAc,SE2$xSEAc,SE3$xSEAc, SE4$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc, SE4$ySEAc)

df_SE <- data.frame(x, y, species)

plot(df_SE$x, df_SE$y)

#plot SE1$xSEAc vs SE1$ySEAc; the ellipse(s) form as a series of points
#This gives SEAc (sample size corrected area)
#hull total area (TA)

find_hull <- function(df) df[chull(df$C, df$N), ]

hulls <- ddply(dat.hair, "species", find_hull)

#Plot data

spps_hulls <- ggplot(dat.hair, aes(x = C, y = N, color = species)) +
  geom_polygon(data = hulls, alpha = 0.1, linetype=3) +
  geom_point(size=3) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  scale_x_continuous(breaks = seq(-30, -20, 1)) +
  scale_y_continuous(breaks = seq(0, 13, 1)) +
  labs(title = "Isotopic Niche", subtitle = "Standard Ellipses and Convex Hulls") +
  geom_path(data=df_SE, aes(x=x, y=y, color=species), linetype=1, size=2) +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = c(0.92, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "18",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "16"))
spps_hulls <- spps_hulls + labs(color = "Mustelid Species")

tiff(file = "figures/hulls_ellipses_spps_18apr2020.tif", 
     units = "cm", width = 30, height = 20, res = 300)

spps_hulls

dev.off()

#calculate convex hulls for mustelids

dat.hairm <- data.matrix(dat.hair, rownames.force = NA)

siber.hull.metrics(dat.hairm$C, dat.hairm$N, dat.hairm$spps_code, R = 10^4)

#calculate % niche overlap between consumers----

data(dat.hair) # 4 species, 2 isotopes

aggregate(dat.hair[8:9], dat.hair[1], mean) # isotope means per species

# format data for plotting function

hair.data <- tapply(1:nrow(dat.hair), dat.hair$species, function(ii) X = dat.hair[ii,8:9])

# niche overlap plots for 95% niche region sizes
# generate parameter draws from the "default" posteriors of each species

nsamples <- 500

system.time({
  hair.par <- tapply(1:nrow(dat.hair), dat.hair$species,
                     function(ii) niw.post(nsamples = nsamples, X = dat.hair[ii,8:9]))
})

# format data for plotting function

hair.data <- tapply(1:nrow(dat.hair), dat.hair$species, function(ii) X = dat.hair[ii,8:9])

clrs <- c("black", "red", "blue", "orange") # colors for each species

niche.plot(niche.par = hair.par, niche.data = hair.data, pfrac = .1,
           iso.names = expression(delta^{15}*N, delta^{13}*C),
           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))

# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.

system.time({
  over <- overlap(hair.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, .99))
})

# posterior expectations of overlap metrics

over.mean <- apply(over*100, c(1:2, 4), mean)

round(over.mean)

write.csv(over.mean, "outputs/over.mean.csv")

# overlap plot

clrs <- c("black", "black", "black", "black") # color for each species

ii <- 1 # which niche region alpha level to use

overlap.plot(over[,,,ii], col = clrs, mean.cred.col = "black",
             xlab = paste0("Overlap Probability (%) -- Niche Region Size: ",
                           dimnames(over)[[4]][ii]))

#plot the SEAc for female mink, male mink, and all marten----
#remove any objects from previous analysis

rm(list = ls())

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

#subset consumer data by species manually

f_mink <- subset(dat.hair, species == 'mink' & sex == 'F')

m_mink <- subset(dat.hair, species == 'mink' & sex == 'M')

marten <- subset(dat.hair, species == 'marten')

otter <- subset(dat.hair, species == 'otter')

#replace mink with f_mink and m_mink in respective dataframes

f_mink <- f_mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "f_mink", as.character(species)))

m_mink <- m_mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "m_mink", as.character(species)))

#make new dataframe with f_mink, m_mink, marten, and otter

mink_marten <- rbind(f_mink, m_mink, marten, otter) 

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(f_mink$C, f_mink$N)
SE2 <- standard.ellipse(m_mink$C, m_mink$N)
SE3 <- standard.ellipse(marten$C, marten$N)
SE4 <- standard.ellipse(otter$C, otter$N)

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

f_mink_ <- rep("f_mink", length(SE1$xSEAc))
m_mink_ <- rep("m_mink", length(SE2$xSEAc))
marten_ <- rep("marten", length(SE3$xSEAc))
otter_ <- rep("otter", length(SE4$xSEAc))

#create new data frame with names and ellipse outputs

species <- c(f_mink_, m_mink_, marten_, otter_)
x <- c(SE1$xSEAc,SE2$xSEAc,SE3$xSEAc, SE4$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc, SE4$ySEAc)

df_SE <- data.frame(x, y, species)

plot(df_SE$x, df_SE$y)

#plot SE1$xSEAc vs SE1$ySEAc; the ellipse(s) form as a series of points
#This gives SEAc (sample size corrected area)
#hull total area (TA)

find_hull <- function(df) df[chull(df$C, df$N), ]

hulls <- ddply(mink_marten, "species", find_hull)

#Plot data

sex_hulls <- ggplot(mink_marten, aes(x = C, y = N, color = species)) +
  geom_polygon(data = hulls, alpha = 0.1, linetype=3) +
  geom_point(size=3) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  scale_x_continuous(breaks = seq(-30, -20, 1)) +
  scale_y_continuous(breaks = seq(0, 13, 1)) +
  labs(title = "Isotopic Niche", subtitle = "Standard Ellipses and Convex Hulls") +
  geom_path(data=df_SE, aes(x=x, y=y, color=species), linetype=1, size=2) +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = c(0.92, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "18",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "16"))
sex_hulls <- sex_hulls + labs(color = "Mustelid Groups")

tiff(file = "figures/hulls_ellipses_sex_18apr2020.tif", 
     units = "cm", width = 30, height = 20, res = 300)

sex_hulls

dev.off()

#calculate convex hulls for mustelids----

siber.hull.metrics(dat.hair$C, dat.hair$N, dat.hair$spps_code, R = 10^4)

#calculate % niche overlap between consumers

data(dat.hair) # 4 species, 2 isotopes

aggregate(dat.hair[8:9], dat.hair[1], mean) # isotope means per species

# format data for plotting function

hair.data <- tapply(1:nrow(dat.hair), dat.hair$species, function(ii) X = dat.hair[ii,8:9])

# niche overlap plots for 95% niche region sizes
# generate parameter draws from the "default" posteriors of each species

nsamples <- 500

system.time({
  hair.par <- tapply(1:nrow(dat.hair), dat.hair$species,
                     function(ii) niw.post(nsamples = nsamples, X = dat.hair[ii,8:9]))
})

# format data for plotting function

hair.data <- tapply(1:nrow(dat.hair), dat.hair$species, function(ii) X = dat.hair[ii,8:9])

clrs <- c("black", "red", "blue", "orange") # colors for each species

niche.plot(niche.par = hair.par, niche.data = hair.data, pfrac = .1,
           iso.names = expression(delta^{15}*N, delta^{13}*C),
           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))

# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.

system.time({
  over <- overlap(hair.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, .99))
})

# posterior expectations of overlap metrics

over.mean <- apply(over*100, c(1:2, 4), mean)

round(over.mean)

write.csv(over.mean, "outputs/over.mean.csv")

# overlap plot

clrs <- c("black", "black", "black", "black") # color for each species

ii <- 1 # which niche region alpha level to use

overlap.plot(over[,,,ii], col = clrs, mean.cred.col = "black",
             xlab = paste0("Overlap Probability (%) -- Niche Region Size: ",
                           dimnames(over)[[4]][ii]))

#calculate % niche overlap between f_mink, m_mink, marten, and otter----

#remove any objects from previous analysis

rm(list = ls())

#load mustelid data if you need to

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

#subset consumer data by species manually

f.mink <- subset(dat.hair, species == 'mink' & sex == 'F')

m.mink <- subset(dat.hair, species == 'mink' & sex == 'M')

marten <- subset(dat.hair, species == 'marten')

otter <- subset(dat.hair, species == 'otter')

#replace mink with f_mink and m_mink in respective dataframes

f.mink <- f.mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "f.mink", as.character(species)))

m.mink <- m.mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "m.mink", as.character(species)))

#make new dataframe with f_mink, m_mink, marten, and otter
dat.sex <- rbind(f.mink, m.mink, marten, otter) 

# isotope means per consumer group

aggregate(dat.sex[8:9], dat.sex[1], mean) 

# format data for plotting function

sex.data <- tapply(1:nrow(dat.sex), dat.sex$species, function(ii) X = dat.sex[ii,8:9])

# niche overlap plots for 95% niche region sizes
# generate parameter draws from the "default" posteriors of each species

nsamples <- 500

system.time({
  sex.par <- tapply(1:nrow(dat.sex), dat.sex$species,
                    function(ii) niw.post(nsamples = nsamples, X = dat.sex[ii,8:9]))
})

# format data for plotting function

dat.sex <- tapply(1:nrow(dat.sex), dat.sex$species, function(ii) X = dat.sex[ii,8:9])

clrs <- c("black", "red", "blue", "orange") # colors for each consumer group

niche.plot(niche.par = sex.par, niche.data = dat.sex, pfrac = .1,
           iso.names = expression(delta^{15}*N, delta^{13}*C),
           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))

# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.

system.time({
  sex <- overlap(sex.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, .99))
})

# posterior expectations of overlap metrics

sex.mean <- apply(sex*100, c(1:2, 4), mean)

round(sex.mean)

write.csv(sex.mean, "outputs/sex.mean.csv")

# overlap plot

clrs <- c("black", "black", "black", "black") # color for each consumer group

ii <- 1 # which niche region alpha level to use

overlap.plot(sex[,,,ii], col = clrs, mean.cred.col = "black",
             xlab = paste0("Overlap Probability (%) -- Niche Region Size: ",
                           dimnames(over)[[4]][ii]))

#Plot prey data as scatter----
#Add prey data and plot as scatter plot
#remove any objects from previous analysis

rm(list = ls())

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot bone from analysis (only muscle)

dat.prey <- dat.prey[-c(29, 30, 31), ]

summary(dat.prey)

#plot prey by species

ggplot(dat.prey, aes(x = C , y = N)) +
  geom_point(aes(shape = species, color = species)) +
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10,11,12,13,
                                14,15,16,17,18,19,20,21,22,23,1,2,3,4,5,6, 7)) +
  scale_color_manual(values = c("violetred", "slateblue4", 
                                "seagreen", "coral3", "chartreuse4", 
                                "chocolate", "chocolate1", "chocolate2",
                                "chocolate3", "coral4", "brown",
                                "darkblue","brown2","brown3","darkmagenta",
                                "deeppink", "darkgoldenrod1", "forestgreen", 
                                "darkgoldenrod3", "gray5", "#FF0000", "#00FF00", 
                                "#3399FF", "#000000", "cyan", "lightskyblue",
                                "indianred", "saddlebrown", "olivedrab",
                                "orange")) +  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Isotopic Signatures", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "16",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.92, 0.98)) +
  scale_x_continuous(breaks = seq(-34, -20, 1)) +
  scale_y_continuous(breaks = seq(0, 13, 1))

#plot prey groups = amphibian, bird, fish, mammal
ggplot(dat.prey, aes(x = C , y = N)) +
  geom_point(aes(shape = group, color = group)) +
  scale_shape_manual(values = c(3,15,16,17)) +
  scale_color_manual(values = c("#FF0000", "#00FF00", "#3399FF", "#000000")) +  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  theme(axis.title.x = element_text(family="Times New Roman", size=11),
        axis.title.y = element_text(family="Times New Roman", size=11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "11"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "11"),
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "11"),
        legend.position = "top") +
  scale_x_continuous(breaks = seq(-34, -20, 1)) +
  scale_y_continuous(breaks = seq(0, 13, 1))

#plot prey group_1 = amphibian, bird, fish, mammal, berry
ggplot(dat.prey, aes(x = C , y = N)) +
  geom_point(aes(shape = group_1, color = group_1)) +
  scale_shape_manual(values = c(3,15,16,17,18)) +
  scale_color_manual(values = c("#FF0000", "#00FF00", "#3399FF", "#000000", "gold1")) +  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  theme(axis.title.x = element_text(family="Times New Roman", size=11),
        axis.title.y = element_text(family="Times New Roman", size=11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "11"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "11"),
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "11"),
        legend.position = "top") +
  scale_x_continuous(breaks = seq(-34, -20, 1)) +
  scale_y_continuous(breaks = seq(0, 13, 1))

#plot mean by prey group_2 = t.vertebrate, fish, berry
#Join means by group_2 to the original data frame and pipe to ggplot

prey_mean <- left_join(dat.prey, 
              dat.prey %>%
               group_by(group_2) %>%
                summarise_at(vars(C, N), funs(mean = mean))
              )%>% 
             ggplot(aes(colour=group_2, shape = group_2)) +

  #Plot raw data
  
  geom_point(aes(x=C, y=N), size = 2, alpha = 0.6) +

  #Plot mean by prey group
  
  geom_point(aes(x=C_mean, y=N_mean), size=4) +

  #Add vertical errorbars 
  
  scale_shape_manual(values = c(15, 21, 17)) +
  scale_color_manual(values = c("gray1", "gray1", "gray9")) +  
  stat_summary(fun.data = mean_se, fun.args = list(mult=1.96), 
               aes(x = C_mean, y=N),
               geom="errorbar", width=0.05) +

  #Add horizontal errorbars
  
  stat_summaryh(fun.data=mean_se_h, fun.args=list(mult=1.96), 
                aes(x=C, y=N_mean),
                geom="errorbarh", width=0.1) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Isotopic Signatures", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "16",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.92, 0.99999)) +
  scale_x_continuous(breaks = seq(-34, -20, 1)) +
  scale_y_continuous(breaks = seq(-5, 13, 1))
  prey_mean <- prey_mean + labs(colour = "Diet Groups", shape = "Diet Groups")

#Export plot

tiff(file = "figures/prey_means_18apr2020.tif", 
     units = "cm", width = 22, height = 17, res = 300)

prey_mean

dev.off()

#K Nearest Neighbor Randomization Test to assign prey groups----

#from: https://www.datacamp.com/community/tutorials/machine-learning-in-r
#acccessed on April 6, 2020

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot bone from analysis (only muscle)

dat.prey <- dat.prey[-c(29, 30, 31), ]

summary(dat.prey)

#Normalize prey N and C values
#Build `normalize()` function

normalize <- function(x) {
  num <- x - min(x)
  denom <- max(x) - min(x)
  return (num/denom)
}

# Normalize the `prey` data

prey.norm <- as.data.frame(lapply(dat.prey[5:6], normalize))

# Summarize `prey.norm`

summary(prey.norm)

#training and test sets to assess model performance later; split data into 2 sets:
#training set (2/3, 0.67) and test set (1/3, 0.33)
#each prey group must have equal chance of being assigned

#set a seed for random number generator

set.seed(1234)

#assign each row 1 (prob. 0.67) or 2 (prob. 0.33), with replacement

ind <- sample(2, nrow(dat.prey), replace=TRUE, prob=c(0.67, 0.33))

# Compose training set

prey.training <- dat.prey[ind==1, 5:6]

# Inspect training set

head(prey.training)

# Compose test set

prey.test <- dat.prey[ind==2, 5:6]

# Inspect test set

head(prey.test)

# Compose `prey` training labels

prey.trainLabels <- dat.prey[ind==1,4]

# Inspect result

print(prey.trainLabels)

# Compose `iris` test labels

prey.testLabels <- dat.prey[ind==2, 4]

# Inspect result

print(prey.testLabels)

# Build the model
#set k = sqrt(n). where n = number of data points in training (should be odd num)
#k = sqrt(N)/2, therefore k = 5

prey_pred <- knn(train = prey.training, test = prey.test, cl = prey.trainLabels, k=5)

# Inspect `prey_pred`

prey_pred

# Put `prey.testLabels` in a data frame

preyTestLabels <- data.frame(prey.testLabels)

# Merge `prey_pred` and `prey.testLabels` 

prey.merge <- data.frame(prey_pred, prey.testLabels)

# Specify column names for `prey.merge`

names(prey.merge) <- c("Predicted Prey Group", "Observed Prey Group")

# Inspect `merge` 

prey.merge

#assess model performance (predicted vs observed)

prey_knn <- CrossTable(x = prey.testLabels, y = prey_pred, prop.chisq=FALSE)

#determine if prey groups are significantly different in C and N----

# Shapiro-Wilk normality test for C13 and N15 by prey group

with(dat.prey, shapiro.test(N[group_2 == "t.vertebrate"]))

with(dat.prey, shapiro.test(C[group_2 == "t.vertebrate"]))

with(dat.prey, shapiro.test(N[group_2 == "berry"]))

with(dat.prey, shapiro.test(C[group_2 == "berry"]))

with(dat.prey, shapiro.test(N[group_2 == "fish"]))

with(dat.prey, shapiro.test(C[group_2 == "fish"]))

#fish N p < 0.05

#kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(C ~ group_2, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

posthoc.kruskal.nemenyi.test(x = C, g = group_2, dist="Tukey")

#t.vertebrates significantly different in C from berry (p = 0.01) and fish (p << 0.001)
#fish not significantly different in C from berry (p = 0.98)

#kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(N ~ group_2, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

posthoc.kruskal.nemenyi.test(x = N, g = group_2, dist="Tukey")

#All groups significantly different in N values (p << 0.001)

#make summary table with prey means and SE----

#Add prey data

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot bone from analysis (only muscle)

dat.prey <- dat.prey[-c(29, 30, 31), ]

#summarize by species(diet item)

summary.prey <- as.data.frame(dat.prey %>%
                                group_by(species) %>%
                                dplyr::summarize(mean_C = mean(C),
                                se_C = sd(C)/sqrt(n()),
                                mean_N = mean(N),
                                se_N=sd(N)/sqrt(n()),
                                sample_size = n()))

summary.prey

write.csv(summary.prey, "outputs/summary_prey_17apr2020.csv")

#Summarize by diet group (group_2)

summary.group <- as.data.frame(dat.prey %>%
                                group_by(group_2) %>%
                                dplyr::summarize(mean_C = mean(C),
                                                 se_C = sd(C)/sqrt(n()),
                                                 mean_N = mean(N),
                                                 se_N=sd(N)/sqrt(n()),
                                                 sample_size = n()))

summary.group

write.csv(summary.group, "outputs/summary_group_17apr2020.csv")

#Summarize consumer stable isotope values (Mean +/- SE) by groupings----

#remove any objects from previous analysis

rm(list = ls())

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

#summarize by consumer species

summary.species <- as.data.frame(dat.hair %>%
                                 group_by(species) %>%
                                 dplyr::summarize(mean_C = mean(C),
                                                  se_C = sd(C)/sqrt(n()),
                                                  mean_N = mean(N),
                                                  se_N=sd(N)/sqrt(n()),
                                                  sample_size = n()))

summary.species

write.csv(summary.species, "outputs/summary_Mustelid_species_17apr2020.csv")

#summarize by f.mink, m.mink, marten, otter

# load SI data (cleaned) from hair of consumers
#spps_code: f.mink = 1, m.mink = 2, marten = 3, and otter = 4

#remove any objects from previous analysis

rm(list = ls())

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

#subset consumer data by species manually

f_mink <- subset(dat.hair, species == 'mink' & sex == 'F')

m_mink <- subset(dat.hair, species == 'mink' & sex == 'M')

marten <- subset(dat.hair, species == 'marten')

otter <- subset(dat.hair, species == 'otter')

#replace mink with f_mink and m_mink in respective dataframes

f_mink <- f_mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "f_mink", as.character(species)))

m_mink <- m_mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "m_mink", as.character(species)))

#fix spps_codes so that f_mink = 1, m_mink = 2, marten = 3, otter = 4

m_mink <- m_mink %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 1, 2, as.character(spps_code)))

marten <- marten %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 2, 3, as.character(spps_code)))

#make new dataframe with f_mink, m_mink, marten, and otter

mink_marten <- rbind(f_mink, m_mink, marten, otter) 

#summarize by f.mink, m.mink, marten, otter

summary.sex <- as.data.frame(mink_marten %>%
                                   group_by(species) %>%
                                   dplyr::summarize(mean_C = mean(C),
                                                    se_C = sd(C)/sqrt(n()),
                                                    mean_N = mean(N),
                                                    se_N=sd(N)/sqrt(n()),
                                                    sample_size = n()))

summary.sex

write.csv(summary.sex, "outputs/summary statistics/summary_sex_17apr2020.csv")

#Run mixing models with mustelid species as grouping variable----

# load SI data (cleaned) from hair of consumers
#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

#remove any objects from previous analysis

rm(list = ls())

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

consumer.spps <- dat.hair[c("spps_code", "C", "N")]

names(consumer.spps)[2] <- "d13C"

names(consumer.spps)[3] <- "d15N"

consumer.matrix <- data.matrix(consumer.spps, rownames.force = NA)

#Add prey data (source data)

dat.prey <- read.csv("data/jprf_prey_data_26mar2020.csv", header = TRUE, sep = ",")

#remove burbot bone from analysis (only muscle)

dat.prey <- dat.prey[-c(29, 30, 31), ]

#make new data frame with Mean and SD for C and N by prey group

dat.prey = data.table(dat.prey)

dat.source = dat.prey[, list("Meand13C" = mean(C),
                       "SDd13C" = sd(C),
                      "Meand15N" = mean(N), 
                      "SDd15N" = sd(N)), 
                   by = c("group_2")]

#convert to matrix; prey groups: t.vertebrate = 3, fish = 2, berry = 1

source.matrix <- data.matrix(dat.source, rownames.force = NA)

#make trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- matrix(c(3, 2, 0.682 , 3, 1.02, 2, 1, 0.341, 2, 0.682, 1, 2, 0.682, 3, 1.02),ncol=5,byrow=TRUE)

colnames(dat.tef) <- c("Source", "Meand13C","SDd13C","Meand15N", "SDd15N")

tef.df <- as.data.frame.matrix(dat.tef)

tef.matrix <- data.matrix(tef.df, rownames.force = NA)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

species.mod <- siarmcmcdirichletv4(consumer.spps, dat.source, tef.df, concdep = 0, 500000, 50000)

#view raw output data

out <- species.mod$output

#summarize model output and write to csv

melted <- melt(out)

summary.species <- ddply(melted, c("Var2"), summarise,
                         mean = mean(value), max = max(value), min = min(value), sd = sd(value),
                         CI = 1.96*(sd(value)/sqrt(length(value))), count = length(value))

write.csv(summary.species, "outputs/summary_species_model.csv")

#plot density distributions by species and diet group----

#density plots for mink

mink_mod <- melted[1:90000,1:3]

# Rename the column and the values in the factor
levels(mink_mod$Var2)[levels(mink_mod$Var2)=="t.vertebrateG1"] <- "T.vertebrate"
levels(mink_mod$Var2)[levels(mink_mod$Var2)=="fishG1"] <- "Fish"
levels(mink_mod$Var2)[levels(mink_mod$Var2)=="berryG1"] <- "Berry"
names(mink_mod)[names(mink_mod)=="Var2"]  <- "Diet_Groups"

mink_diet <- ggplot(data=mink_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,.8), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  ggtitle("Mink Diet Estimates") +
  labs(x = "Proportion of Diet", y = "Density") +
        theme(title = element_text(family="Times New Roman", size=16),
              axis.title.x = element_text(family="Times New Roman", size=16),
              axis.title.y = element_text(family="Times New Roman", size=16),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
              axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
              legend.text = element_text(family = "Times New Roman", colour = "black", size=16),
              legend.position = "right")
              
tiff(file = "figures/density distributions/by species/mink_diet_18apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

mink_diet

dev.off()

#for marten

marten_mod <- melted[150001:240000,1:3]

levels(marten_mod$Var2)[levels(marten_mod$Var2)=="t.vertebrateG2"] <- "T.vertebrate"
levels(marten_mod$Var2)[levels(marten_mod$Var2)=="fishG2"] <- "Fish"
levels(marten_mod$Var2)[levels(marten_mod$Var2)=="berryG2"] <- "Berry"
names(marten_mod)[names(marten_mod)=="Var2"]  <- "Diet_Groups"

marten_diet <- ggplot(data=marten_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,11), expand = c(0, 0)) +
  ggtitle("Marten Diet Estimates") +
  labs(x = "Proportion of Diet", y = "Density") +
  theme(title = element_text(family="Times New Roman", size=16),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = "right")

tiff(file = "figures/density distributions/by species/marten_diet_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

marten_diet

dev.off()

#for weasel

weasel_mod <- melted[300001:390000,1:3]

# Rename the column and the values in the factor
levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="t.vertebrateG3"] <- "T.vertebrate"
levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="fishG3"] <- "Fish"
levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="berryG3"] <- "Berry"
names(weasel_mod)[names(weasel_mod)=="Var2"]  <- "Diet_Groups"

weasel_diet <- ggplot(data=weasel_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,15), expand = c(0, 0)) +
  ggtitle("Weasel Diet Estimates") +
  labs(x = "Proportion of Diet", y = "Density") +
  theme(title = element_text(family="Times New Roman", size=16),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = "right")

tiff(file = "figures/density distributions/by species/weasel_diet_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

weasel_diet

dev.off()

#for otter

otter_mod <- melted[450001:540000,1:3]

# Rename the column and the values in the factor

levels(otter_mod$Var2)[levels(otter_mod$Var2)=="t.vertebrateG4"] <- "T.vertebrate"
levels(otter_mod$Var2)[levels(otter_mod$Var2)=="fishG4"] <- "Fish"
levels(otter_mod$Var2)[levels(otter_mod$Var2)=="berryG4"] <- "Berry"
names(otter_mod)[names(otter_mod)=="Var2"]  <- "Diet_Groups"

otter_diet <- ggplot(data=otter_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  ggtitle("Otter Diet Estimates") +
  labs(x = "Proportion of Diet", y = "Density") +
  theme(title = element_text(family="Times New Roman", size=16),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.8, 0.99))

tiff(file = "figures/density distributions/by species/otter_diet_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

otter_diet

dev.off()

#combine into one figure

mink_diet <- ggplot(data=mink_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "Diet Estimates", subtitle = "Mink") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

marten_diet <- ggplot(data=marten_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,15), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "", subtitle = "Marten") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = "none")

weasel_diet <- ggplot(data=weasel_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,15), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "", subtitle = "Weasel") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

otter_diet <- ggplot(data=otter_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "", subtitle = "Otter") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.78, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 16,face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
        otter_diet <- otter_diet + guides(fill=guide_legend(title="Diet Groups"))

        otter_diet        
tiff("figures/density distributions/by species/all_species_18apr2020.tiff", units="cm", width=30, height=20, res=300)

ggarrange(mink_diet, otter_diet, marten_diet, weasel_diet, 
          labels = c("A", "B", "C", "D"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=18, face = "bold", color = "black"),
          ncol = 2, nrow = 2)

dev.off()

#plot density distributions by diet item with overlapping consumer spps----

melted$Var3 <- ifelse(grepl("t.vertebrate", melted$Var2), "t.vertebrate", "NA")

mink <- subset(melted, grepl("G1", melted$Var2))

mink$species  <- "mink"

marten <- subset(melted, grepl("G2", melted$Var2))

marten$species  <- "marten"

weasel <- subset(melted, grepl("G3", melted$Var2))

weasel$species  <- "weasel"

otter <- subset(melted, grepl("G4", melted$Var2))

otter$species  <- "otter"

#make new dataframe with all species

spps_density <- rbind(mink, marten, weasel, otter) 

#plot density distributions by species for t.vertebrate

spps_t.vert <- subset(spps_density, Var3 == 't.vertebrate')

t.vert_spps <- ggplot(data=spps_t.vert, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,6), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Terrestrial vertebrate") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.91, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
t.vert_spps <- t.vert_spps + guides(fill=guide_legend(title="Mustelid Species"))

tiff(file = "figures/density distributions/by species/spps by diet groups/t.vert_spps_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

t.vert_spps

dev.off()

#plot density distributions by species for fish

spps_density$Var4 <- ifelse(grepl("fish", melted$Var2), "fish", "NA")

spps_fish <- subset(spps_density, Var4 == 'fish')

fish_spps <- ggplot(data=spps_fish, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,8), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Fish") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.91, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
fish_spps <- fish_spps + guides(fill=guide_legend(title="Mustelid Species"))

tiff(file = "figures/density distributions/by species/spps by diet groups/fish_spps_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

fish_spps

dev.off()

#plot density distributions by species for fish

spps_density$Var4 <- ifelse(grepl("berry", melted$Var2), "berry", "NA")

spps_berry <- subset(spps_density, Var4 == 'berry')

berry_spps <- ggplot(data=spps_berry, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Berry") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.91, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
berry_spps <- berry_spps + guides(fill=guide_legend(title="Mustelid Species"))

tiff(file = "figures/density distributions/by species/spps by diet groups/berry_spps_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

berry_spps

dev.off()

#combine probability densities into one plot
#with 3 columns and 1 row

t.vert_spps <- ggplot(data=spps_t.vert, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "Diet Group", subtitle = "Terrestrial vertebrate") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

fish_spps <- ggplot(data=spps_fish, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Fish") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

berry_spps <- ggplot(data=spps_berry, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Berry") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.73, 0.95),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
berry_spps <- berry_spps + guides(fill=guide_legend(title="Mustelid Species"))

tiff("figures/density distributions/by species/spps by diet groups/all_prey_17apr2020.tiff", units="cm", width=40, height=10, res=300)

ggarrange(t.vert_spps, fish_spps, berry_spps,
          labels = c("A", "B", "C"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=16, face = "bold", color = "black"),
          ncol = 3, nrow = 1)

dev.off()

#arrange in 1 column with 3 rows

t.vert_spps <- ggplot(data=spps_t.vert, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "Diet Group", subtitle = "Terrestrial vertebrate") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5,1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.72, 0.95),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "18",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
t.vert_spps <- t.vert_spps + guides(fill=guide_legend(title="Mustelid Species"))

fish_spps <- ggplot(data=spps_fish, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Fish") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

berry_spps <- ggplot(data=spps_berry, aes(x= value, group = species, fill = species)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Berry") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

tiff("figures/density distributions/by species/spps by diet groups/all_prey_1col_18apr2020.tiff", units="cm", width=15, height=30, res=300)

ggarrange(t.vert_spps, fish_spps, berry_spps,
          labels = c("A", "B", "C"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=18, face = "bold", color = "black"),
          ncol = 1, nrow = 3)

dev.off()

#Run mixing models with f.mink, m.mink, marten, and otter as consumer groups----

# load SI data (cleaned) from hair of consumers
#spps_code: f.mink = 1, m.mink = 2, marten = 3, and otter = 4

#remove any objects from previous analysis

rm(list = ls())

dat.hair <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.hair$C, plot = FALSE)$out

dat.hair[which(dat.hair$C %in% outlier),]

dat.hair <- dat.hair[-which(dat.hair$C %in% outlier),]

#subset consumer data by species manually

f_mink <- subset(dat.hair, species == 'mink' & sex == 'F')

m_mink <- subset(dat.hair, species == 'mink' & sex == 'M')

marten <- subset(dat.hair, species == 'marten')

otter <- subset(dat.hair, species == 'otter')

#replace mink with f_mink and m_mink in respective dataframes

f_mink <- f_mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "f_mink", as.character(species)))

m_mink <- m_mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "m_mink", as.character(species)))

#fix spps_codes so that f_mink = 1, m_mink = 2, marten = 3, otter = 4

m_mink <- m_mink %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 1, 2, as.character(spps_code)))

marten <- marten %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 2, 3, as.character(spps_code)))

#make new dataframe with f_mink, m_mink, marten, and otter

mink_marten <- rbind(f_mink, m_mink, marten, otter) 

#format consumer data as matrix for SIAR

consumer.sex <- mink_marten[c("spps_code", "C", "N")]

names(consumer.sex)[2] <- "d13C"

names(consumer.sex)[3] <- "d15N"

consumer.matrix <- data.matrix(consumer.sex, rownames.force = NA)

#Add prey data (source data)

dat.prey <- read.csv("data/jprf_prey_data_26mar2020.csv", header = TRUE, sep = ",")

#remove burbot bone from analysis (only muscle)

dat.prey <- dat.prey[-c(29, 30, 31), ]

#make new data frame with Mean and SD for C and N by prey group

dat.prey = data.table(dat.prey)

dat.source = dat.prey[, list("Meand13C" = mean(C),
                             "SDd13C" = sd(C),
                             "Meand15N" = mean(N), 
                             "SDd15N" = sd(N)), 
                      by = c("group_2")]

#convert to matrix; prey groups: t.vertebrate = 3, fish = 2, berry = 1

source.matrix <- data.matrix(dat.source, rownames.force = NA)

#make trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- matrix(c(3, 2, 0.682 , 3, 1.02, 2, 1, 0.341, 2, 0.682, 1, 2, 0.682, 3, 1.02),ncol=5,byrow=TRUE)

colnames(dat.tef) <- c("Source", "Meand13C","SDd13C","Meand15N", "SDd15N")

tef.df <- as.data.frame.matrix(dat.tef)

tef.matrix <- data.matrix(tef.df, rownames.force = NA)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

sex.mod <- siarmcmcdirichletv4(consumer.sex, dat.source, tef.df, concdep = 0, 500000, 50000)

#save model output and write to csv

sex.out <- sex.mod$output

sex.melted <- melt(sex.out)

summary.sex <- ddply(sex.melted, c("Var2"), summarise,
                         mean = mean(value), max = max(value), min = min(value), sd = sd(value),
                         CI = 1.96*(sd(value)/sqrt(length(value))), count = length(value))

sex.is.num <- sapply(summary.sex, is.numeric)
summary.sex[sex.is.num] <- lapply(summary.sex[sex.is.num], round, 2)


summary.sex


write.csv(summary.sex, "outputs/mixing model outputs/summary_sex_model_19apr2020.csv")

#plot density distributions by sex (mink) and diet group----

#density plots for female mink

f.mink_mod <- sex.melted[1:90000,1:3]

levels(f.mink_mod$Var2)[levels(f.mink_mod$Var2)=="t.vertebrateG1"] <- "T.vertebrate"
levels(f.mink_mod$Var2)[levels(f.mink_mod$Var2)=="fishG1"] <- "Fish"
levels(f.mink_mod$Var2)[levels(f.mink_mod$Var2)=="berryG1"] <- "Berry"
names(f.mink_mod)[names(f.mink_mod)=="Var2"]  <- "Diet_Groups"

f.mink_diet <- ggplot(data=f.mink_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "Mustelid Group", subtitle = "Female mink") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.9, 0.99999),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
f.mink_diet <- f.mink_diet + guides(fill=guide_legend(title="Diet Groups"))

tiff(file = "figures/density distributions/by sex (mink)/f.mink_diet_17apr2020.tif", 
      units = "cm", width = 20, height = 12, res = 300)

f.mink_diet

dev.off()

#for male mink

m.mink_mod <- sex.melted[150001:240000,1:3]

levels(m.mink_mod$Var2)[levels(m.mink_mod$Var2)=="t.vertebrateG2"] <- "T.vertebrate"
levels(m.mink_mod$Var2)[levels(m.mink_mod$Var2)=="fishG2"] <- "Fish"
levels(m.mink_mod$Var2)[levels(m.mink_mod$Var2)=="berryG2"] <- "Berry"
names(m.mink_mod)[names(m.mink_mod)=="Var2"]  <- "Diet_Groups"

m.mink_diet <- ggplot(data=m.mink_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "Mustelid Group", subtitle = "Male mink") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.9, 0.99999),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
m.mink_diet <- m.mink_diet + guides(fill=guide_legend(title="Diet Groups"))

tiff(file = "figures/density distributions/by sex (mink)/m.mink_diet_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

m.mink_diet

dev.off()

#for marten

marten_mod <- sex.melted[300001:390000,1:3]

levels(marten_mod$Var2)[levels(marten_mod$Var2)=="t.vertebrateG3"] <- "T.vertebrate"
levels(marten_mod$Var2)[levels(marten_mod$Var2)=="fishG3"] <- "Fish"
levels(marten_mod$Var2)[levels(marten_mod$Var2)=="berryG3"] <- "Berry"
names(marten_mod)[names(marten_mod)=="Var2"]  <- "Diet_Groups"

marten_diet <- ggplot(data=marten_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "Mustelid Group", subtitle = "Marten") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.9, 0.99999),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
marten_diet <- marten_diet + guides(fill=guide_legend(title="Diet Groups"))

tiff(file = "figures/density distributions/by sex (mink)/marten_diet_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

marten_diet

dev.off()

#for otter

otter_mod <- sex.melted[450001:540000,1:3]

levels(otter_mod$Var2)[levels(otter_mod$Var2)=="t.vertebrateG4"] <- "T.vertebrate"
levels(otter_mod$Var2)[levels(otter_mod$Var2)=="fishG4"] <- "Fish"
levels(otter_mod$Var2)[levels(otter_mod$Var2)=="berryG4"] <- "Berry"
names(otter_mod)[names(otter_mod)=="Var2"]  <- "Diet_Groups"

otter_diet <- ggplot(data=otter_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "Mustelid Group", subtitle = "Otter") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.9, 0.99999),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
otter_diet <- otter_diet + guides(fill=guide_legend(title="Diet Groups"))

tiff(file = "figures/density distributions/by sex (mink)/otter_diet_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

otter_diet

dev.off()

#combine densities into one plot by Mustelid group (f.mink, m.mink, marten, otter)

f.mink_diet <- ggplot(data=f.mink_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "Mustelid Group", subtitle = "Female mink") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

m.mink_diet <- ggplot(data=m.mink_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "", subtitle = "Male mink") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.76, 0.99999),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 18, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
m.mink_diet <- m.mink_diet + guides(fill=guide_legend(title="Diet Groups"))


marten_diet <- ggplot(data=marten_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "", subtitle = "Marten") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

otter_diet <- ggplot(data=otter_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density", title = "", subtitle = "Otter") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

tiff("figures/density distributions/by sex (mink)/all_sex_18apr2020.tiff", units="cm", width=30, height=20, res=300)

ggarrange(f.mink_diet, m.mink_diet, marten_diet, otter_diet,
          labels = c("A", "B", "C", "D"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=18, face = "bold", color = "black"),
          ncol = 2, nrow = 2)

dev.off()

#plot density distributions by diet item with overlapping consumer spps

sex.melted$Var3 <- ifelse(grepl("t.vertebrate", sex.melted$Var2), "t.vertebrate", "NA")

f.mink <- subset(sex.melted, grepl("G1", sex.melted$Var2))

f.mink$group  <- "F.mink"

m.mink <- subset(sex.melted, grepl("G2", sex.melted$Var2))

m.mink$group  <- "M.mink"

marten <- subset(sex.melted, grepl("G3", sex.melted$Var2))

marten$group  <- "Marten"

otter <- subset(sex.melted, grepl("G4", sex.melted$Var2))

otter$group  <- "Otter"

#make new dataframe with all species

sex_density <- rbind(f.mink, m.mink, marten, otter) 

#plot density distributions by species for t.vertebrate

sex_t.vert <- as.data.frame(subset(sex_density, Var3 == 't.vertebrate'))

t.vert_sex <- ggplot(data=sex_t.vert, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,7), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "Diet Group", subtitle = "Terrestrial vertebrate") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.91, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
t.vert_sex <- t.vert_sex + guides(fill=guide_legend(title="Mustelid Groups"))

tiff(file = "figures/density distributions/by sex (mink)/sex by diet group/t.vert_sex_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

t.vert_sex

dev.off()

#plot density distributions by species for fish

sex_density$Var4 <- ifelse(grepl("fish", sex.melted$Var2), "fish", "NA")

sex_fish <- subset(sex_density, Var4 == 'fish')

fish_sex <- ggplot(data=sex_fish, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,8), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Fish") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.91, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))

fish_sex <- fish_sex + guides(fill=guide_legend(title="Mustelid Groups"))

tiff(file = "figures/density distributions/by sex (mink)/sex by diet group/fish_sex_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

fish_sex

dev.off()

#plot density distributions by species for fish

sex_density$Var4 <- ifelse(grepl("berry", sex.melted$Var2), "berry", "NA")

sex_berry <- subset(sex_density, Var4 == 'berry')

berry_sex <- ggplot(data=sex_berry, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Berry") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.91, 0.99),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
berry_sex <- berry_sex + guides(fill=guide_legend(title="Mustelid Groups"))

tiff(file = "figures/density distributions/by sex (mink)/sex by diet group/berry_sex_17apr2020.tif", 
     units = "cm", width = 20, height = 12, res = 300)

berry_sex

dev.off()

#combine probability densities into one plot
#with 3 columns and 1 row

t.vert_sex <- ggplot(data=sex_t.vert, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "Diet Group", subtitle = "Terrestrial vertebrate") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

fish_sex <- ggplot(data=sex_fish, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Fish") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

berry_sex <- ggplot(data=sex_berry, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Berry") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.73, 0.95),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "14",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
berry_sex <- berry_sex + guides(fill=guide_legend(title="Mustelid Groups"))

tiff("figures/density distributions/by sex (mink)/sex by diet group/all_prey_sex_17apr2020.tiff", units="cm", width=40, height=10, res=300)

ggarrange(t.vert_sex, fish_sex, berry_sex,
          labels = c("A", "B", "C"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=16, face = "bold", color = "black"),
          ncol = 3, nrow = 1)

dev.off()

#arrange in 1 column with 3 rows

t.vert_sex <- ggplot(data=sex_t.vert, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "Diet Group", subtitle = "Terrestrial vertebrate") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position = c(0.72, 0.95),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 18, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=16))
t.vert_sex <- t.vert_sex + guides(fill=guide_legend(title="Mustelid Groups"))

fish_sex <- ggplot(data=sex_fish, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Fish") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

berry_sex <- ggplot(data=sex_berry, aes(x= value, group = group, fill = group)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  labs(x = "Proportion of Diet", y = "Density", title = "", subtitle = "Berry") +
  theme(plot.title = element_text(family="Times New Roman", size=18, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=17, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=16),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=16),
        legend.position="none",
        legend.text = element_blank())

tiff("figures/density distributions/by sex (mink)/sex by diet group/all_prey__sex_1col_18apr2020.tiff", 
          units="cm", width=15, height=30, res=300)

ggarrange(t.vert_sex, fish_sex, berry_sex,
          labels = c("A", "B", "C"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=18, face = "bold", color = "black"),
          ncol = 1, nrow = 3)

dev.off()

#compare nail and hair SI values among marten and mink----
#load SI data (cleaned) from hair and nails of consumers
#spps_code: mink = 1

dat.mink <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                         sql = "select * from file where `species` = 'mink'", sep = ",")
#make tissue code


#plot mean and 95% CI for mink by tissue type

minktis_means <- left_join(dat.mink, 
                            dat.mink %>%
                              group_by(Tissue) %>%
                              summarise_at(vars(C, N), funs(mean = mean))
) %>% 
  ggplot(aes(colour=Tissue, shape = Tissue)) +
  # Plot raw data
  geom_point(aes(x=C, y=N), size = 2, alpha = 0.7) +
  # Plot mean by species
  geom_point(aes(x=C_mean, y=N_mean), size=3) +
  # Add vertical errorbars 
  scale_shape_manual(values = c(0, 21, 17)) +
  scale_color_manual(values = c("gray1", "gray1", "gray15")) + 
  stat_summary(fun.data = mean_se, fun.args = list(mult=1.96), 
               aes(x = C_mean, y=N),
               geom="errorbar", width=0.05) +
  # Add horizontal errorbars
  stat_summaryh(fun.data=mean_se_h, fun.args=list(mult=1.96), 
                aes(x=C, y=N_mean),
                geom="errorbarh", width=0.1) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Isotopic Signatures", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "16",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.88, 0.98)) +
  scale_x_continuous(breaks = seq(-30, -20, 1)) +
  scale_y_continuous(breaks = seq(0, 14, 1))
minktis_means <- minktis_means + labs(colour = "Mink Tissues", shape = "Mink Tissues")

tiff(file = "figures/mink_tissue_95CI_28apr2020.tif", 
     units = "cm", width = 20, height = 15, res = 300)

minktis_means

dev.off()

#do the same for marten tissues (spps_code: marten = 2)

dat.marten <- read.csv.sql("data/jprf_consumer_data_27mar2020.csv", header = TRUE, 
                           sql = "select * from file where `species` = 'marten'", sep = ",")


#Find and remove outliers and blank data from dataset (NA = 9999)

outlier <- boxplot(dat.marten$C, plot = FALSE)$out

dat.marten[which(dat.marten$C %in% outlier),]

dat.marten <- dat.marten[-which(dat.marten$C %in% outlier),]

outlier <- boxplot(dat.marten$N, plot = FALSE)$out

dat.marten[which(dat.marten$N %in% outlier),]

dat.marten <- dat.marten[-which(dat.marten$N %in% outlier),]

boxplot(dat.marten$C)

boxplot(dat.marten$N)

#plot mean and 95% CI for marten by tissue type

marttis_means <- left_join(dat.marten, 
                           dat.marten %>%
                             group_by(Tissue) %>%
                             summarise_at(vars(C, N), funs(mean = mean))
) %>% 
  ggplot(aes(colour=Tissue, shape = Tissue)) +
  # Plot raw data
  geom_point(aes(x=C, y=N), size = 2, alpha = 0.7) +
  # Plot mean by species
  geom_point(aes(x=C_mean, y=N_mean), size=3) +
  # Add vertical errorbars 
  scale_shape_manual(values = c(0, 21, 17)) +
  scale_color_manual(values = c("gray1", "gray1", "gray15")) + 
  stat_summary(fun.data = mean_se, fun.args = list(mult=1.96), 
               aes(x = C_mean, y=N),
               geom="errorbar", width=0.05) +
  # Add horizontal errorbars
  stat_summaryh(fun.data=mean_se_h, fun.args=list(mult=1.96), 
                aes(x=C, y=N_mean),
                geom="errorbarh", width=0.1) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Isotopic Signatures", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "16",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.88, 0.98)) +
  scale_x_continuous(breaks = seq(-30, -20, 1)) +
  scale_y_continuous(breaks = seq(0, 14, 1))
marttis_means <- marttis_means + labs(colour = "Marten Tissues", shape = "Marten Tissues")

tiff(file = "figures/marten_tissue_95CI_28apr2020.tif", 
     units = "cm", width = 20, height = 15, res = 300)

marttis_means

dev.off()

# Shapiro-Wilk normality test for C13 and N15 by species

with(dat.mink, shapiro.test(N[Tissue == "hair"]))

with(dat.mink, shapiro.test(C[Tissue == "hair"]))

with(dat.mink, shapiro.test(N[Tissue == "nailb"]))

with(dat.mink, shapiro.test(C[Tissue == "nailb"]))

with(dat.mink, shapiro.test(N[Tissue == "nailt"]))

with(dat.mink, shapiro.test(C[Tissue == "nailt"]))


###mink nailt (tip of nail) SI data do not follow a normal distribution (p < 0.01)

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(C ~ Tissue, data = dat.mink)

require(PMCMR)

data(dat.mink)

attach(dat.mink)

posthoc.kruskal.nemenyi.test(x = C, g = "Tissue", dist="Tukey")

###no diff in C btwn mink and marten (p = .1), but diff btwn mink & weasel (p << 0.001)
###       and btwn marten and weasel (p = .02) and btwn weasel and otter (p = 0.004)
##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(N ~ spps_code, data = dat.hair)

require(PMCMR)

data(dat.hair)

attach(dat.hair)

posthoc.kruskal.nemenyi.test(x = N, g = spps_code, dist="Tukey")

