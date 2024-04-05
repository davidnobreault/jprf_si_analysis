##These are the SI analyses of diet of mink, marten, weasel, and otter in JPRF ----
#performed by D. Breault in 2019-2020 for Dexter - last updated 10jun2020

#set library

.libPaths("C:/R/library")

#set working directory

setwd(("C:/R/projects/jprf_si_analysis"))

#remove any objects from previous analysis

rm(list = ls())

options("install.lock" = FALSE)

#Install packages ----

install.packages("devtools")

install.packages("extrafont")

install.packages("ggpubr")

install.packages("PMCMRplus")

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

install.packages("ggrepel")

install.packages("descr")

install.packages("car")

install.packages('splancs')

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

library('PMCMRplus')

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

library("ggrepel")

library("descr")

library("car")

library("splancs")

#TABLE 1. Summarize consumer stable isotope values (Mean +/- SE) by species, and tissue----

#remove any objects from previous analysis

rm(list = ls())

#load SI data (cleaned) from consumers

dat.mustel <- read.csv("data/jprf_consumer_data_7jun2021.csv", header = TRUE, sep = ",")

#make values in Tissue column either hair or nail

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and marten

dat.mustel <- setNames(aggregate(list(dat.mustel$C, dat.mustel$N), 
                                by = list(dat.mustel$ID, dat.mustel$Tissue, dat.mustel$spps_code, 
                                          dat.mustel$species, dat.mustel$sex), mean),
                              c("ID", "Tissue", "spps_code", "species", "sex", "C", "N"))

#subset consumer data by species manually

mink.hair <- subset(dat.mustel, species == 'mink' & Tissue == 'hair')

mink.nail <- subset(dat.mustel, species == 'mink' & Tissue == 'nail')



marten.hair <- subset(dat.mustel, species == 'marten' & Tissue == 'hair')

marten.nail <- subset(dat.mustel, species == 'marten' & Tissue == 'nail')

weasel <- subset(dat.mustel, species == 'weasel')

otter <- subset(dat.mustel, species == 'otter')

#replace species with different unique codes in respective dataframes

marten.hair <- data.frame(append(marten.hair, c(tissue='marten.hair')))

marten.nail <- data.frame(append(marten.nail, c(tissue='marten.nail')))

mink.hair <- data.frame(append(mink.hair, c(tissue='mink.hair')))

mink.nail <- data.frame(append(mink.nail, c(tissue='mink.nail')))

otter <- data.frame(append(otter, c(tissue='otter')))

weasel <- data.frame(append(weasel, c(tissue='weasel')))

#make new dataframe with all mustelid groups for comparison

dat.mustel.all <- rbind(marten.hair, marten.nail, 
                        mink.hair, mink.nail, weasel, otter) 

#summarize by mink_n, and mink_h

summary.mustel.all <- as.data.frame(dat.mustel.all %>%
                               group_by(tissue) %>%
                               dplyr::summarize(sample_size = n(),
                                                mean_C = mean(C),
                                                se_C = sd(C)/sqrt(n()),
                                                range_C = range (C),
                                                mean_N = mean(N),
                                                se_N=sd(N)/sqrt(n()),
                                                range_N = range (N)))

summary.mustel.all

write.csv(summary.mustel.all, "outputs/summary statistics/summary_mustelids_all_13dec2022.csv")

#TABLE 1. Niche Ellipses with mustelid species and tissues as consumer groupings----
#remove any objects from previous analysis

rm(list = ls())

#load SI data (cleaned) from consumers

dat.mustel <- read.csv("data/jprf_consumer_data_7jun2021.csv", header = TRUE, sep = ",")

#make values in Tissue column either hair or nail

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and marten

dat.means <- setNames(aggregate(list(dat.mustel$C, dat.mustel$N), 
                                by = list(dat.mustel$ID, dat.mustel$Tissue, dat.mustel$spps_code, 
                                          dat.mustel$species, dat.mustel$sex), mean),
                      c("ID", "Tissue", "spps_code", "species", "sex", "C", "N"))

#subset consumer data by species manually

mink_h <- subset(dat.means, species == 'mink' & Tissue == 'hair')

mink_n <- subset(dat.means, species == 'mink' & Tissue == 'nail')

marten_h <- subset(dat.means, species == 'marten'& Tissue == 'hair')

marten_n <- subset(dat.means, species == 'marten'& Tissue == 'nail')

weasel <- subset(dat.means, species == 'weasel')

otter <- subset(dat.means, species == 'otter')

#replace mink with mink_h and mink_n in respective dataframes

mink_h <- mink_h %>% 
  mutate(species = ifelse(as.character(species) == "mink", "mink_h", as.character(species)))

mink_n <- mink_n %>% 
  mutate(species = ifelse(as.character(species) == "mink", "mink_n", as.character(species)))

marten_h <- marten_h %>% 
  mutate(species = ifelse(as.character(species) == "marten", "marten_h", as.character(species)))

marten_n <- marten_n %>% 
  mutate(species = ifelse(as.character(species) == "marten", "marten_n", as.character(species)))

#use mustelid data to calculate standard ellipse metrics (SEAc) by species

minkh_SEA <- data.frame(standard.ellipse(mink_h$C, mink_h$N, confs = NULL, steps = 1))

write.csv(minkh_SEA, "outputs/standard.ellipse statistics/minkh_SEA_08jun2021.csv")

minkn_SEA <- data.frame(standard.ellipse(mink_n$C, mink_n$N, confs = NULL, steps = 1))

write.csv(minkn_SEA, "outputs/standard.ellipse statistics/minkn_SEA_08jun2021.csv")

martenh_SEA <- data.frame(standard.ellipse(marten_h$C, marten_h$N, confs = NULL, steps = 1))

write.csv(martenh_SEA, "outputs/standard.ellipse statistics/martenh_SEA_08jun2021.csv")

martenn_SEA <- data.frame(standard.ellipse(marten_n$C, marten_n$N, confs = NULL, steps = 1))

write.csv(martenn_SEA, "outputs/standard.ellipse statistics/martenn_SEA_08jun2021.csv")

weasel_SEA <- data.frame(standard.ellipse(weasel$C, weasel$N, confs = NULL, steps = 2))

write.csv(weasel_SEA, "outputs/standard.ellipse statistics/weasel_SEA_08jun2021.csv")

otter_SEA <- standard.ellipse(otter$C, otter$N, confs = NULL, steps = 1)

write.csv(otter_SEA, "outputs/standard.ellipse statistics/otter_SEA_08jun2021.csv")

#Otter niche metrics:
#SEA = 4.345662, SEAc = 4.562946, theta = -0.447343

#Unpaired, two sample t-test: niche breadth of larger consumer greater than small consumer----

#remove any objects from previous analysis

rm(list = ls())

#load data

dat.body <- read.csv("data/body_size.csv", header = TRUE, sep = ",")

#plot

ggscatter(dat.body, x = "Weight", y = "SEAc", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Weight (g)", ylab = "Niche Breadth")

#Spearman rank correlation

body.rank <-cor.test(dat.body$Weight, dat.body$SEAc,  method = "spearman")

body.rank

#are mink and marten significantly diff weights? If so, which one is heavier?

#remove any objects from previous analysis

rm(list = ls())

#load data

dat.weight <- read.csv("data/weight_peltoff_adult.csv", header = TRUE, sep = ",")

#subset consumer data by species manually

marten_mink <- dat.weight[1:40,]

marten_mink

#shapiro-wilks normality test

with(marten_mink, shapiro.test(weight_g[ï..species == "Marten"]))

with(marten_mink, shapiro.test(weight_g[ï..species == "Mink"]))

### no violations in normality (p > 0.05)

# Levene's test with one independent variable to test for homogeneity of variances

levtestweight <- leveneTest(data = marten_mink, y = marten_mink$weight_g, group = marten_mink$ï..species)

levtestweight

#no violation of homegeneity of variances

# Compute t-test

res <- t.test(weight_g ~ ï..species, data = marten_mink, var.equal = TRUE)
res

#plot weight vs niche breadth (SEAc)

#Figure 1 - weight (g) versus niche breadths

dat.weight

weight_niche <- as.data.frame(dat.weight %>%
                             group_by(ï..species) %>%
                             dplyr::summarize(mean_w = mean(weight_g),
                                              se_w = sd(weight_g)/sqrt(n()),
                                              mean_n = mean(sea)))

bodysize_niche <- ggplot(data = weight_niche, aes(colour=ï..species, shape = ï..species)) +

  # Plot mean by species
  geom_point(aes(x=mean_w, y=mean_n), size=2) +

  # Add horizontal errorbars
  stat_summary(fun.data = mean_se_h, fun.args=list(mult=1.96), 
                aes(x=mean_w, y=mean_n),
                geom="errorbar", width=0.1) +
  ylab("Dietary Niche Breadth (SEAc)") +
  xlab("Weight (g)") +
  theme(plot.title = element_text(family="Times New Roman", size=9, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=8, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=9),
        axis.title.y = element_text(family="Times New Roman", size=9),
        panel.background = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "8"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "8"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "9",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "8"),
        legend.key = element_blank(),
        legend.position = c(0.91, 0.93),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_continuous(breaks = seq(0, 10000, 1000), limits = c(0,10000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0,10), expand = c(0, 0))

bodysize_niche <- bodysize_niche + labs(fill = "Mustelid Species", 
                                        colour = "Mustelid Species", shape = "Mustelid Species")

tiff(file = "figures/figure2_31July2023.tif", 
     units = "mm", width = 174, height = 139, res = 300)

bodysize_niche

dev.off()

#Are mustelid species (& mink m and f) significantly different in C & N?----

#remove any objects from previous analysis

rm(list = ls())

#load SI data (cleaned) from consumers

dat.mustel <- read.csv("data/jprf_consumer_data_7jun2021.csv", header = TRUE, sep = ",")

#Take mean of bottom and top nail segments for mink and marten

dat.means <- setNames(aggregate(list(dat.mustel$C, dat.mustel$N), 
                                by = list(dat.mustel$ID, dat.mustel$Tissue, dat.mustel$spps_code, 
                                          dat.mustel$species), mean),
                      c("ID", "Tissue", "spps_code", "species", "C", "N"))

#subset consumer data by species manually

mink_h <- subset(dat.means, species == 'mink' & Tissue == 'hair')

marten_h <- subset(dat.means, species == 'marten'& Tissue == 'hair')

weasel <- subset(dat.means, species == 'weasel')

otter <- subset(dat.means, species == 'otter')

#replace mink with f_mink and m_mink in respective dataframes

mink_h <- mink_h %>% 
  mutate(species = ifelse(as.character(species) == "mink", "mink_h", as.character(species)))

marten_h <- marten_h %>% 
  mutate(species = ifelse(as.character(species) == "marten", "marten_h", as.character(species)))

#make new dataframe with mink_h, marten_h, weasel, and otter

dat.mustel.all <- rbind(mink_h, marten_h, weasel, otter) 

# Shapiro-Wilk normality test for C13 and N15 by species

with(dat.mustel.all, shapiro.test(N[species == "mink_h"]))

with(dat.mustel.all, shapiro.test(C[species == "mink_h"]))

with(dat.mustel.all, shapiro.test(N[species == "marten_h"]))

with(dat.mustel.all, shapiro.test(C[species == "marten_h"]))

with(dat.mustel.all, shapiro.test(N[species == "weasel"]))

with(dat.mustel.all, shapiro.test(C[species == "weasel"]))

with(dat.mustel.all, shapiro.test(N[species == "otter"]))

with(dat.mustel.all, shapiro.test(C[species == "otter"]))

### no violations in normality (p > 0.01)

# Levene's test with one independent variable to test for homogeneity of variances

levtestC <- leveneTest(data = dat.mustel.all, y = dat.mustel.all$C, group = dat.mustel.all$species)

levtestC

levtestN <- leveneTest(data = dat.mustel.all, y = dat.mustel.all$N, group = dat.mustel.all$species)

levtestN

#violations in assumption of equal variances (p = 0.002)

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(C ~ species, data = dat.mustel.all)

dat.mustel.all$species <- as.factor(dat.mustel.all$species)

kwAllPairsNemenyiTest( C ~ species, data = dat.mustel.all)

#data: C by species

#       marten_h mink_h otter  
#mink_h 0.0316   -      -      
#otter  8.3e-07  0.0086 -      
#weasel 0.9976   0.0272 9.5e-07

pairwise.wilcox.test(dat.mustel.all$C, dat.mustel.all$species,
                     p.adjust.method = "BH")

#       marten_h mink_h  otter  
#mink_h 0.00603  -       -      
#otter  1.6e-06  0.00041 -      
#weasel 0.89339  0.00338 1.7e-07

kruskal.test(N ~ species, data = dat.mustel.all)

kwAllPairsNemenyiTest( N ~ species, data = dat.mustel.all)

#data: N by species

#       marten_h mink_h otter  
#mink_h 0.0136   -      -      
#otter  2.5e-08  0.0022 -      
#weasel 0.9785   0.0694 5.9e-07

pairwise.wilcox.test(dat.mustel.all$N, dat.mustel.all$species,
                     p.adjust.method = "BH")

#       marten_h mink_h  otter  
#mink_h 0.0023   -       -      
#otter  7.4e-11  5.2e-05 -      
#weasel 0.6010   0.0094  2.9e-09

#is N of otter greater than N of mink?

wilcox.test(mink_h$N, otter$N, alternative="less")

#is N of marten greater than N of weasel?

wilcox.test(marten_h$N, weasel$N)


#FIGURE 1. Plot the carnivore means and 95% CI, by species, and tissue (mink/marten)----

#remove any objects from previous analysis

rm(list = ls())

#load SI data (cleaned) from consumers

dat.mustel <- read.csv("data/jprf_consumer_data_7jun2021.csv", header = TRUE, sep = ",")

#make values in Tissue column either hair or nail

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and marten

dat.means <- setNames(aggregate(list(dat.mustel$C, dat.mustel$N), 
                                by = list(dat.mustel$ID, dat.mustel$Tissue, dat.mustel$spps_code, 
                                          dat.mustel$species), mean),
                      c("ID", "Tissue", "spps_code", "species", "C", "N"))

#subset consumer data by species manually

mink_h <- subset(dat.means, species == 'mink' & Tissue == 'hair')

mink_n <- subset(dat.means, species == 'mink' & Tissue == 'nail')

marten <- subset(dat.means, species == 'marten'& Tissue == 'hair')

weasel <- subset(dat.means, species == 'weasel')

otter <- subset(dat.means, species == 'otter')

#make new dataframe with all mustelid hair data

dat.mustel.hair <- rbind(mink_h, marten, weasel, otter) 

#Figure 1A. Join means by tissues for mustelids to the original data frame and pipe to ggplot

consumer_means <- left_join(dat.mustel.hair, 
                            dat.mustel.hair %>%
                              group_by(species) %>%
                              summarise_at(vars(C, N), funs(mean = mean))
) %>% 
  ggplot(aes(colour=species, shape = species)) +
  # Plot raw data
  geom_point(aes(x=C, y=N), size = 1.7, alpha = 0.8) +
  # Plot mean by species
  geom_point(aes(x=C_mean, y=N_mean), size=2) +
  # Add vertical errorbars 
  scale_shape_manual(values = c(15, 17, 18, 19)) +
  scale_color_manual(values = c("red", "blue", "green", "orange")) + 
  stat_summary(fun.data = mean_se, fun.args = list(mult=1.96), 
               aes(x = C_mean, y=N),
               geom="errorbar", width=0.05) +
  # Add horizontal errorbars
  stat_summaryh(fun.data=mean_se_h, fun.args=list(mult=1.96), 
                aes(x=C, y=N_mean),
                geom="errorbarh", width=0.1) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  #labs(colour = "Mustelid Groups", shape = "Mustelid Groups") +
  #labs(title = "Isotopic Signatures", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=9, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=8, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=9),
        axis.title.y = element_text(family="Times New Roman", size=9),
        panel.background = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "8"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "8"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "9",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "8"),
        legend.key = element_blank(),
        legend.position = c(0.91, 0.93),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_continuous(breaks = seq(-32, -20, 1)) + 
  scale_y_continuous(breaks = seq(6, 14, 1))

consumer_means <- consumer_means + labs(fill = "Mustelid Species", 
                                        colour = "Mustelid Species", shape = "Mustelid Species")

tiff(file = "figures/figure1_29jun2022.tif", 
     units = "mm", width = 174, height = 139, res = 300)

consumer_means

dev.off()

#Are nail and hair SI values different among marten and mink?----
#load SI data (cleaned) from hair and nails of consumers
#spps_code: mink = 1

#remove any objects from previous analysis

rm(list = ls())

dat.mink <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `species` = 'mink'", sep = ",")

#make values in Tissue column either hair or nail

dat.mink <- dat.mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mink <- dat.mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

# Shapiro-Wilk normality test for C13 and N15 by tissue type among mink

with(dat.mink, shapiro.test(N[Tissue == "hair"]))

with(dat.mink, shapiro.test(C[Tissue == "hair"]))

with(dat.mink, shapiro.test(N[Tissue == "nail"]))

with(dat.mink, shapiro.test(C[Tissue == "nail"]))

with(dat.mink, shapiro.test(N[Tissue == "nailt"]))

with(dat.mink, shapiro.test(C[Tissue == "nailt"]))

with(dat.mink, shapiro.test(N[Tissue == "nailb"]))

with(dat.mink, shapiro.test(C[Tissue == "nailb"]))

# Levene's test with one independent variable to test for homogeneity of variances

levtestC <- leveneTest(data = dat.mink, y = dat.mink$C, group = dat.mink$Tissue)

levtestC

levtestN <- leveneTest(data = dat.mink, y = dat.mink$N, group = dat.mink$Tissue)

levtestN

levtestC <- leveneTest(data = dat.mink, y = dat.mink$C, group = dat.mink$Tissue)

levtestC

levtestN <- leveneTest(data = dat.mink, y = dat.mink$N, group = dat.mink$Tissue)

levtestN

# Compute the analysis of variance

minkC.aov <- aov(C ~ Tissue, data = dat.mink)

minkN.aov <- aov(N ~ Tissue, data = dat.mink)

# Summary of the analysis

summary(minkC.aov)

summary(minkN.aov)

#summary(minkC.aov)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Tissue       1   1.67   1.666   0.454  0.504
#Residuals   54 198.37   3.674            
 
#summary(minkN.aov)
#            Df Sum Sq Mean Sq F value Pr(>F)
#Tissue       1   1.61   1.608   0.724  0.398
#Residuals   54 119.91   2.220  

#remove any objects from previous analysis

rm(list = ls())

dat.marten <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `species` = 'marten'", sep = ",")

#make values in Tissue column either hair or nail

dat.marten <- dat.marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.marten <- dat.marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

# Shapiro-Wilk normality test for C13 and N15 from mink tissues

with(dat.marten, shapiro.test(N[Tissue == "hair"]))

with(dat.marten, shapiro.test(C[Tissue == "hair"]))

with(dat.marten, shapiro.test(N[Tissue == "nail"]))

with(dat.marten, shapiro.test(C[Tissue == "nail"]))

# Levene's test with one independent variable to test for homogeneity of variances

levtestC <- leveneTest(y = dat.marten$C, group = dat.marten$Tissue, data = dat.marten)

levtestC

levtestN <- leveneTest(y = dat.marten$N, group = dat.marten$Tissue, data = dat.marten)

levtestN

#Mann-Whitney-Wilcoxon Test of differences between 2 groups in C13

wilcox.test(C ~ Tissue, data = dat.marten) 

#Mann-Whitney-Wilcoxon Test of differences between 2 groups in N15

wilcox.test(N ~ Tissue, data = dat.marten) 

#kruskal-wallis test of differences between 2 groups in C13

dat.marten$Tissue <- as.factor(dat.marten$Tissue)

kruskal.test(C ~ Tissue, data = dat.marten)

#data:  C by Tissue
#Kruskal-Wallis chi-squared = 9.525, df = 1, p-value = 0.002027

kruskal.test(N ~ Tissue, data = dat.marten)

#data:  N by Tissue
#Kruskal-Wallis chi-squared = 40.365, df = 1, p-value = 2.107e-10

#FIGURE 2. Plot ALL CARNIVORE HAIR ellipses and convex hulls by species with prey (mean +/- SE)----

rm(list = ls())

dat.hair <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

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

SE1
SE2
SE3
SE4

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

mink_ <- rep("mink", length(SE1$xSEAc))
marten_ <- rep("marten", length(SE2$xSEAc))
weasel_ <- rep("weasel", length(SE3$xSEAc))
otter_ <- rep("otter", length(SE4$xSEAc))

#create new data frame with names and ellipse outputs

species <- c(mink_, marten_, weasel_, otter_)
x <- c(SE1$xSEAc,SE2$xSEAc, SE3$xSEAc, SE4$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc, SE3$ySEAc, SE4$ySEAc)

df_SE <- data.frame(x, y, species)

plot(df_SE$x, df_SE$y)

#plot SE1$xSEAc vs SE1$ySEAc; the ellipse(s) form as a series of points
#This gives SEAc (sample size corrected area)
#hull total area (TA)

find_hull <- function(df) df[chull(df$C, df$N), ]

hulls <- ddply(dat.hair, "species", find_hull)

#Add prey data and plot 

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, salamanders, toads, lake trout, salmon, kokanee, squirrel, grouse from analysis

dat.prey <- dat.prey[-c(4:6, 22:23, 26:31, 44:49, 57:60, 67:70), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Squirrel", "squirrel", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Clam", "mussel", as.character(group_2)))

#correct prey by TEF values as in mixing models
#dC13 = +2 for terrestrial vertebrates and berries
#       +1 for fish, clams
#dN15 = +3 for terrestrial vertebrates and berries
#       +2 for fish, clams

dat.prey$C <- dat.prey$C + ifelse(dat.prey$group_2 == "fish.waterfowl", 1, 
                                  ifelse(dat.prey$group_2 == "mussel", 1, 2))

dat.prey$N <- dat.prey$N + ifelse(dat.prey$group_2 == "fish.waterfowl", 2, 
                                  ifelse(dat.prey$group_2 == "mussel", 2, 3))

#Combine the mustelid ellipses with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(group_2) %>%
                             dplyr::summarize(mean_C = mean(C),
                                              se_C = sd(C)/sqrt(n()),
                                              mean_N = mean(N),
                                              se_N=sd(N)/sqrt(n()),
                                              sample_size = n()))

mustelid_w_prey <- ggplot(NULL) +
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = group_2), size = 2) +
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_path(data = df_SE, aes(x = x, y = y, color = species), linetype = 1, size = 1) +
  
  scale_shape_manual(values = c(1,2,0,5,6,7)) +
  scale_color_manual(values = c("red", "blue", "green", "orange")) + 
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  #labs(title = "Mustelid Hair with Prey", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=9, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=8, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=9),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "8"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "8"),
        legend.position = c(0.9, .225), legend.direction = "vertical",
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "9",face = "bold"),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "8"),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")) +
          guides(color = guide_legend(order=1),
          shape = guide_legend(order=2)) +
  scale_x_continuous(breaks = seq(-33, -18, 1)) +
  scale_y_continuous(breaks = seq(-0, 13, 1))

mustelid_w_prey <- mustelid_w_prey + labs(shape = "Prey", color = "Mustelid Species")

tiff(file = "figures/figure2_29jun2022.tif", 
     units = "mm", width = 174, height = 139, res = 300)

mustelid_w_prey

dev.off()
                           
#Figure 3. Plot Mink and Marten Hair and Nail ellipses with prey (Mean +/- SE)----

rm(list = ls())

dat.mima <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `species` = 'mink' 
                         or `species` = 'marten'", sep = ",")

#make values in Tissue column either hair or nail

dat.mima <- dat.mima %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mima <- dat.mima %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and marten

dat.mima <- setNames(aggregate(list(dat.mima$C, dat.mima$N), 
                                 by = list(dat.mima$ID, dat.mima$Tissue, dat.mima$spps_code, 
                                           dat.mima$species, dat.mima$sex), mean),
                       c("ID", "Tissue", "spps_code", "species", "sex", "C", "N"))

#subset consumer data by species (mink and marten) manually

mink <- subset(dat.mima, species == 'mink')

marten <- subset(dat.mima, species == 'marten')

#subset consumer data by tissue type manually

mink.hair <- subset(mink, Tissue == 'hair')

mink.nail <- subset(mink, Tissue == 'nail')

marten.hair <- subset(marten, Tissue == 'hair')

marten.nail <- subset(marten, Tissue == 'nail')

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(mink.hair$C, mink.hair$N)
SE2 <- standard.ellipse(mink.nail$C, mink.nail$N)
SE3 <- standard.ellipse(marten.hair$C, marten.hair$N)
SE4 <- standard.ellipse(marten.nail$C, marten.nail$N)

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

mink.hair_ <- rep("mink.hair", length(SE1$xSEAc))
mink.nail_ <- rep("mink.nail", length(SE2$xSEAc))
marten.hair_ <- rep("marten.hair", length(SE3$xSEAc))
marten.nail_ <- rep("marten.nail", length(SE4$xSEAc))

#create new data frame with names and ellipse outputs

tissues <- c(mink.hair_, mink.nail_, marten.hair_, marten.nail_)
x <- c(SE1$xSEAc,SE2$xSEAc, SE3$xSEAc, SE4$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc,SE3$ySEAc)

df_SE <- data.frame(x, y, tissues)

plot(df_SE$x, df_SE$y)

#Add prey data and plot 

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, clam, salamanders, toads, lake trout, salmon, kokanee, grouse from analysis

dat.prey <- dat.prey[-c(4:6, 22:23, 26:31, 35:38, 44:49, 57:60, 67:70), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Squirrel", "squirrel", as.character(group_2)))

#correct prey by TEF values as in mixing models
#dC13 = +2 for terrestrial vertebrates and berries
#       +1 for fish
#dN15 = +3 for terrestrial vertebrates and berries
#       +2 for fish

dat.prey$C <- dat.prey$C + ifelse(dat.prey$group_2 == "fish.waterfowl", 1, 2)

dat.prey$N <- dat.prey$N + ifelse(dat.prey$group_2 == "fish.waterfowl", 2, 3)

#Combine the mink ellipses and hulls with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(group_2) %>%
                             dplyr::summarize(mean_C = mean(C),
                                              se_C = sd(C)/sqrt(n()),
                                              mean_N = mean(N),
                                              se_N=sd(N)/sqrt(n()),
                                              sample_size = n()))

mima_tissues_w_prey <- ggplot(NULL) +
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = group_2), size = 2) +


  geom_path(data = df_SE, aes(x = x, y = y, color = tissues), linetype = 1, size = 1) +
  
  scale_shape_manual(values = c(1,0,5,6,7)) +
  scale_color_manual(values = c("blue", "red", "green", "orange")) + 
  
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  #labs(title = "Mink and Marten Tissues with Prey", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=9, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=8, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=9),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "8"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "8"),
        legend.position = c(0.87, .2), legend.direction = "vertical",
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "9",face = "bold"),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "8"),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2)) +
  scale_x_continuous(breaks = seq(-33, -18, 1)) +
  scale_y_continuous(breaks = seq(-0, 13, 1))

mima_tissues_w_prey <- mima_tissues_w_prey + labs(shape = "Prey", color = "Mink and Marten Tissues")

tiff(file = "figures/figure3_27july2021.tif", 
     units = "mm", width = 174, height = 139, res = 300)

mima_tissues_w_prey

dev.off()

#TABLE 2. calculate % niche overlap between consumers incl. marten nail and hair----

rm(list = ls())

#load SI data (cleaned) from consumers

dat.mustel <- read.csv("data/jprf_consumer_data_7jun2021.csv", header = TRUE, sep = ",")

#make values in Tissue column either hair or nail

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and marten

dat.mustel <- setNames(aggregate(list(dat.mustel$C, dat.mustel$N), 
                                by = list(dat.mustel$ID, dat.mustel$Tissue, dat.mustel$spps_code, 
                                          dat.mustel$species, dat.mustel$sex), mean),
                      c("ID", "Tissue", "spps_code", "species", "sex", "C", "N"))

#subset consumer data by species manually

marten.hair <- subset(dat.mustel, species == 'marten' & Tissue == 'hair')

marten.nail <- subset(dat.mustel, species == 'marten' & Tissue == 'nail')

mink.hair <- subset(dat.mustel, species == 'mink' & Tissue == 'hair')

mink.nail <- subset(dat.mustel, species == 'mink' & Tissue == 'nail')

weasel <- subset(dat.mustel, species == 'weasel')

otter <- subset(dat.mustel, species == 'otter')

#replace species with different unique codes in respective dataframes

marten.hair <- data.frame(append(marten.hair, c(tissue_spps='marten.hair')))

marten.nail <- data.frame(append(marten.nail, c(tissue_spps='marten.nail')))

mink.hair <- data.frame(append(mink.hair, c(tissue_spps='mink.hair')))

mink.nail <- data.frame(append(mink.nail, c(tissue_spps='mink.nail')))

weasel <- data.frame(append(weasel, c(tissue_spps='weasel')))

otter <- data.frame(append(otter, c(tissue_spps='otter')))

#make new dataframe with all mustelid groups for comparison

dat.mustel.all <- rbind(marten.hair, marten.nail,
                        mink.hair, mink.nail, otter, weasel) 

# format data for plotting function

dat.mustel.all$tissue_spps <- as.factor(dat.mustel.all$tissue_spps)

# niche overlap plots for 95% niche region sizes
# generate parameter draws from the "default" posteriors of each species

nsamples <- 500

system.time({
  mustel.par <- tapply(1:nrow(dat.mustel.all), dat.mustel.all$tissue_spps,
                     function(ii) niw.post(nsamples = nsamples, X = dat.mustel.all[ii,6:7]))
})

# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.

system.time({
  over <- overlap(mustel.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, .99))
})

# posterior expectations of overlap metrics

over.mean <- apply(over*100, c(1:2, 4), mean)

round(over.mean)

write.csv(over.mean, "outputs/over.mean_tissue_spps_28nov2022.csv")

over.95ci <- apply(over*100, c(1:2, 4), )

#Overlap plot

clrs <- c("black", "red", "blue", "orange", "green", "purple") # colors for each species
over.stat <- overlap(mustel.par, nreps = nsamples, nprob = 1e3, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

#TABLE 3. make summary table with prey species means and SE----

rm(list = ls())

#Add prey data

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove salmon, kokanee, beaver, grouse, salamanders, toads, from analysis

dat.prey <- dat.prey[-c(4:6, 22:23, 29:31, 44:46, 57:60, 67:70), ]

#summarize by species(diet item)

summary.prey <- as.data.frame(dat.prey %>%
                                group_by(species) %>%
                                dplyr::summarize(mean_C = mean(C),
                                                 se_C = sd(C)/sqrt(n()),
                                                 mean_N = mean(N),
                                                 se_N=sd(N)/sqrt(n()),
                                                 sample_size = n()))

summary.prey

write.csv(summary.prey, "outputs/summary_prey_18Jun2020.csv")

#K Nearest Neighbor Randomization Test to assign MINK prey groups----

#from: https://www.datacamp.com/community/tutorials/machine-learning-in-r
#acccessed on April 6, 2020

#Add prey data and plot as scatter plot

rm(list = ls())

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, clam, salamanders, toads, lake trout, salmon, kokanee, squirrel, grouse from analysis

dat.prey <- dat.prey[-c(4:6, 16:18, 22:23, 26:31, 35:38, 44:49, 57:60, 67:70), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

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
#k = sqrt(N)/2, therefore k = 3

prey_pred <- knn(train = prey.training, test = prey.test, cl = prey.trainLabels, k=3)

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

prey_knn <- CrossTable(x = prey.testLabels, y = prey_pred, prop.chisq = FALSE)

print(prey_knn)

#determine if MINK prey groups are significantly different in C and N----

# Shapiro-Wilk normality test for C13 and N15 by prey group

with(dat.prey, shapiro.test(N[group_2 == "smallmamm.bird"]))

with(dat.prey, shapiro.test(C[group_2 == "smallmamm.bird"]))

with(dat.prey, shapiro.test(N[group_2 == "fish.waterfowl"]))

with(dat.prey, shapiro.test(C[group_2 == "fish.waterfowl"]))

with(dat.prey, shapiro.test(N[group_2 == "berry"]))

with(dat.prey, shapiro.test(C[group_2 == "berry"]))

with(dat.prey, shapiro.test(N[group_2 == "muskrat"]))

with(dat.prey, shapiro.test(C[group_2 == "muskrat"]))

#kruskal-wallis test of differences between > 2 groups in C13

dat.prey$group_2 <- as.factor(dat.prey$group_2)

kruskal.test(C ~ group_2, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

kwAllPairsNemenyiTest(data = dat.prey, x = C, g = group_2, dist="Tukey")

#               berry fish.waterfowl muskrat
#fish.waterfowl 0.999 -              -      
#muskrat        0.031 0.029          -      
#smallmamm.bird 0.037 0.013          0.406  

#kruskal-wallis test of differences between > 2 groups in N15

dat.prey$group_2 <- as.factor(dat.prey$group_2)

kruskal.test(N ~ group_2, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

kwAllPairsNemenyiTest(data = dat.prey, x = N, g = group_2, dist="Tukey")

#               berry   fish.waterfowl muskrat
#fish.waterfowl 1.1e-10 -              -      
#muskrat        0.37548 0.45616        -      
#smallmamm.bird 0.00874 0.00025        0.99862      

#Run mixing models with female and male mink (hair) as grouping variable----

# load SI data (cleaned) from hair of consumers
#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

#remove any objects from previous analysis

rm(list = ls())

dat.hair <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#subset consumer data by species manually

f.mink <- subset(dat.hair, species == 'mink' & sex == 'F')

m.mink <- subset(dat.hair, species == 'mink' & sex == 'M')

#replace mink with f_mink and m_mink in respective dataframes

f_mink <- f.mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "f_mink", as.character(species)))

m_mink <- m.mink %>% 
  mutate(species = ifelse(as.character(species) == "mink", "m_mink", as.character(species)))

#fix spps_codes so that f_mink = 1, m_mink = 2, marten = 3, weasel = 4, otter = 5

m_mink <- m_mink %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 1, 2, as.character(spps_code)))

#make new dataframe with f_mink, m_mink, marten, and otter

dat.mink.sex <- rbind(f_mink, m_mink) 

consumer.spps <- dat.mink.sex[c("spps_code", "C", "N")]

names(consumer.spps)[2] <- "d13C"

names(consumer.spps)[3] <- "d15N"

consumer.matrix <- data.matrix(consumer.spps, rownames.force = NA)

#Add prey data (source data)

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, clam, salamanders, toads, lake trout, salmon, kokanee, squirrel, grouse from analysis

dat.prey <- dat.prey[-c(4:6, 16:18, 22:23, 26:31, 35:38, 44:49, 57:60, 67:70), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

#make new data frame with Mean and SD for C and N by prey group
#convert group_2 to numeric: berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "berry", 1, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "fish.waterfowl", 2, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "smallmamm.bird", 3, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "muskrat", 4, as.character(group_2)))

dat.prey = data.table(dat.prey)

dat.source = dat.prey[, list("Meand13C" = mean(C),
                       "SDd13C" = sd(C),
                      "Meand15N" = mean(N), 
                      "SDd15N" = sd(N)), 
                   by = c("group_2")]

#convert to matrix

source.matrix <- data.matrix(dat.source)

#make trophic enrichment fractor table with TEF values for 4 sources in this model
#berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4
#dC13 = +2 for terrestrial vertebrates and berries
#       +1 for fish, clams
#dN15 = +3 for terrestrial vertebrates and berries
#       +2 for fish, clams

dat.tef <- matrix(c(1, 2, 0.682, 3, 1.02, 
                    2, 1, 0.341, 2, 0.682, 
                    3, 2, 0.682, 3, 1.02,
                    4, 2, 0.682, 3, 1.02), ncol=5,byrow=TRUE)

colnames(dat.tef) <- c("Source", "Meand13C","SDd13C","Meand15N", "SDd15N")

tef.df <- as.data.frame.matrix(dat.tef)

tef.matrix <- data.matrix(tef.df, rownames.force = NA)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

species.mod <- siarmcmcdirichletv4(consumer.matrix, source.matrix, tef.matrix, concdep = 0, 500000, 50000)

#view raw output data

out <- species.mod$output

melted <- melt(out)

#TABLE. summarize model output and write to csv----

summary.species <- ddply(melted, c("Var2"), summarise,
                         mean = mean(value), max = max(value), min = min(value), sd = sd(value),
                         CI = 1.96*(sd(value)/sqrt(length(value))), count = length(value))

write.csv(summary.species, "outputs/mixing model outputs/summary_mink_sex_08jun2021.csv")

#ESM 1. plot density distributions by female and male mink----

#density plots for female mink

fmink_mod <- melted[1:120000,1:3]

mmink_mod <- melted[180001:300000,1:3]

# Save an object to a file

save(species.mod, melted, fmink_mod, mmink_mod, file = "RData/sex_mink_9jun2021.RData")

# Restore the objects

load(file = "RData/sex_mink_9jun2021.RData")

# Rename the column and the values in the factor
#convert group_2 to numeric: #berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4

levels(fmink_mod$Var2)[levels(fmink_mod$Var2)=="1G1"] <- "berry"
levels(fmink_mod$Var2)[levels(fmink_mod$Var2)=="2G1"] <- "fish.waterfowl"
levels(fmink_mod$Var2)[levels(fmink_mod$Var2)=="3G1"] <- "smallmamm.bird"
levels(fmink_mod$Var2)[levels(fmink_mod$Var2)=="4G1"] <- "muskrat"
names(fmink_mod)[names(fmink_mod)=="Var2"]  <- "Diet_Groups"

levels(mmink_mod$Var2)[levels(mmink_mod$Var2)=="1G2"] <- "berry"
levels(mmink_mod$Var2)[levels(mmink_mod$Var2)=="2G2"] <- "fish.waterfowl"
levels(mmink_mod$Var2)[levels(mmink_mod$Var2)=="3G2"] <- "smallmamm.bird"
levels(mmink_mod$Var2)[levels(mmink_mod$Var2)=="4G2"] <- "muskrat"
names(mmink_mod)[names(mmink_mod)=="Var2"]  <- "Diet_Groups"

fmink_diet <- ggplot(data=fmink_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=9, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=8),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=8),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=8))
fmink_diet <- fmink_diet + guides(fill=guide_legend(title="Diet Groups"))

mmink_diet <- ggplot(data=mmink_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(plot.title = element_text(family="Times New Roman", size=9, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=8, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=9, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=8),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=8),
        legend.position="none",
        legend.text = element_blank())

tiff(file = "figures/supplementary material/ESM_1_27july2021.tif", 
     units = "mm", width = 156, height = 234, res = 300)

ggarrange(fmink_diet, mmink_diet, 
          labels = c("A", "B"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=10, face = "bold", color = "black"),
          ncol = 1, nrow = 2)

dev.off()

#Do female mink eat more terrestrial prey than male? Do male eat more aquatic?----

# Restore the objects

load(file = "RData/sex_grizzly_long2.RData")

#make new columns with consumer groups and diet item names

fmink_mod <- data.frame(append(fmink_mod, c(group ='fmink')))

mmink_mod <- data.frame(append(mmink_mod, c(group ='mmink')))

levels(fmink_mod$Var2)[levels(fmink_mod$Var2)=="1G1"] <- "berry"
levels(fmink_mod$Var2)[levels(fmink_mod$Var2)=="2G1"] <- "fish.waterfowl"
levels(fmink_mod$Var2)[levels(fmink_mod$Var2)=="3G1"] <- "smallmamm.bird"
levels(fmink_mod$Var2)[levels(fmink_mod$Var2)=="4G1"] <- "muskrat"
names(fmink_mod)[names(fmink_mod)=="Var2"]  <- "Diet_Groups"

levels(mmink_mod$Var2)[levels(mmink_mod$Var2)=="1G2"] <- "berry"
levels(mmink_mod$Var2)[levels(mmink_mod$Var2)=="2G2"] <- "fish.waterfowl"
levels(mmink_mod$Var2)[levels(mmink_mod$Var2)=="3G2"] <- "smallmamm.bird"
levels(mmink_mod$Var2)[levels(mmink_mod$Var2)=="4G2"] <- "muskrat"
names(mmink_mod)[names(mmink_mod)=="Var2"]  <- "Diet_Groups"

#make new dataframe with all mustelid groups for comparison

dat.all <- rbind(fmink_mod, mmink_mod) 


p.meat.F.M <- 

p.meat.F.M <- as.data.frame(p.meat.F.M)

#Run mixing models with mink hair and nail as grouping variable----

# load SI data (cleaned) from hair of consumers
#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

#remove any objects from previous analysis

rm(list = ls())

dat.mink <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                           sql = "select * from file where `species` = 'mink'", sep = ",")

#make values in Tissue column either hair or nail

dat.mink <- dat.mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mink <- dat.mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and mink

dat.means <- setNames(aggregate(list(dat.mink$C, dat.mink$N), 
                                by = list(dat.mink$ID, dat.mink$Tissue, dat.mink$spps_code, 
                                          dat.mink$species, dat.mink$sex), mean),
                      c("ID", "Tissue", "spps_code", "species", "sex", "C", "N"))

#subset consumer data by tissue type manually

mink_h <- subset(dat.means, Tissue == 'hair')

mink_n <- subset(dat.means, Tissue == 'nail')

#replace labels in respective dataframes

mink_h <- mink_h %>% 
  mutate(species = ifelse(as.character(species) == "mink", "mink_h", as.character(species)))

mink_n <- mink_n %>% 
  mutate(species = ifelse(as.character(species) == "mink", "mink_n", as.character(species)))

#fix spps_codes so that mink_h = 1, mink_n = 2

mink_n <- mink_n %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 1, 2, as.character(spps_code)))

#make new dataframe 

dat.mink.tis <- rbind(mink_h, mink_n) 

consumer.spps <- dat.mink.tis[c("spps_code", "C", "N")]

names(consumer.spps)[2] <- "d13C"

names(consumer.spps)[3] <- "d15N"

consumer.matrix <- data.matrix(consumer.spps, rownames.force = NA)

#Add prey data (source data)

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, clam, salamanders, toads, lake trout, salmon, kokanee, squirrel, grouse from analysis

dat.prey <- dat.prey[-c(4:6, 16:18, 22:23, 26:31, 35:38, 44:49, 57:60, 67:70), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

#make new data frame with Mean and SD for C and N by prey group
#convert group_2 to numeric: berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "berry", 1, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "fish.waterfowl", 2, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "smallmamm.bird", 3, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "muskrat", 4, as.character(group_2)))

dat.prey = data.table(dat.prey)

dat.source = dat.prey[, list("Meand13C" = mean(C),
                             "SDd13C" = sd(C),
                             "Meand15N" = mean(N), 
                             "SDd15N" = sd(N)), 
                      by = c("group_2")]

#convert to matrix

source.matrix <- data.matrix(dat.source)

#make trophic enrichment fractor table with TEF values for 4 sources in this model
#berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4
#dC13 = +2 for berry, smallmamm.bird, squirrel, and muskrat
#       +1 for fish.waterfowl, clam
#dN15 = +3 for berry, smallmamm.bird, squirrel, and muskrat
#       +2 for fish.waterfowl, clam

dat.tef <- matrix(c(1, 2, 0.682, 3, 1.02, 
                    2, 1, 0.341, 2, 0.682, 
                    3, 2, 0.682, 3, 1.02,
                    4, 2, 0.682, 3, 1.02), ncol=5,byrow=TRUE)

colnames(dat.tef) <- c("Source", "Meand13C","SDd13C","Meand15N", "SDd15N")

tef.df <- as.data.frame.matrix(dat.tef)

tef.matrix <- data.matrix(tef.df, rownames.force = NA)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

species.mod <- siarmcmcdirichletv4(consumer.matrix, source.matrix, tef.matrix, concdep = 0, 500000, 50000)

#view raw output data

out <- species.mod$output

#TABLE. summarize model output and write to csv----

melted <- melt(out)

summary.species <- ddply(melted, c("Var2"), summarise,
                         mean = mean(value), max = max(value), min = min(value), sd = sd(value),
                         CI = 1.96*(sd(value)/sqrt(length(value))), 
                         quant_5 = quantile(value, 0.05), 
                         quant_95 = quantile(value, 0.95), 
                         count = length(value))

write.csv(summary.species, "outputs/mixing model outputs/summary_mink_tissues_10nov2022.csv")

#ESM 2. plot density distributions by mink tissues----

#density plots for female mink

minkh_mod <- melted[1:120000,1:3]

minkn_mod <- melted[180001:300000,1:3]

# Save an object to a file

save(species.mod, melted, minkh_mod, minkn_mod, file = "RData/tissue_mink_27july2021.RData")

# Restore the objects

load(file = "RData/tissue_mink_27july2021.RData")

# Rename the column and the values in the factor
#convert group_2 to numeric: #berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4

levels(minkh_mod$Var2)[levels(minkh_mod$Var2)=="1G1"] <- "berry"
levels(minkh_mod$Var2)[levels(minkh_mod$Var2)=="2G1"] <- "fish.waterfowl"
levels(minkh_mod$Var2)[levels(minkh_mod$Var2)=="3G1"] <- "smallmamm.bird"
levels(minkh_mod$Var2)[levels(minkh_mod$Var2)=="4G1"] <- "muskrat"
names(minkh_mod)[names(minkh_mod)=="Var2"]  <- "Diet_Groups"

levels(minkn_mod$Var2)[levels(minkn_mod$Var2)=="1G2"] <- "berry"
levels(minkn_mod$Var2)[levels(minkn_mod$Var2)=="2G2"] <- "fish.waterfowl"
levels(minkn_mod$Var2)[levels(minkn_mod$Var2)=="3G2"] <- "smallmamm.bird"
levels(minkn_mod$Var2)[levels(minkn_mod$Var2)=="4G2"] <- "muskrat"
names(minkn_mod)[names(minkn_mod)=="Var2"]  <- "Diet_Groups"

minkh_diet <- ggplot(data=minkh_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,40), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=9, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=8),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=8),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=8))
minkh_diet <- minkh_diet + guides(fill=guide_legend(title="Diet Groups"))

minkn_diet <- ggplot(data=minkn_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=9, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=8),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=8),
        legend.position="none",
        legend.text = element_blank())

tiff(file = "figures/supplementary material/ESM_2_27july2021.tif", 
     units = "mm", width = 156, height = 234, res = 300)

ggarrange(minkh_diet, minkn_diet, 
          labels = c("A", "B"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=10, face = "bold", color = "black"),
          ncol = 1, nrow = 2)

dev.off()

#K Nearest Neighbor Randomization Test to assign MARTEN/WEASEL prey groups----

#from: https://www.datacamp.com/community/tutorials/machine-learning-in-r
#acccessed on April 6, 2020

#Add prey data and plot as scatter plot

rm(list = ls())

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove fish, beaver, muskrat, grouse, waterfowl, salamanders, and toads from analysis

dat.prey <- dat.prey[-c(4:6, 22:74), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Squirrel", "squirrel", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

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
#k = sqrt(N)/2, therefore k = 3

prey_pred <- knn(train = prey.training, test = prey.test, cl = prey.trainLabels, k=3)

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

prey_knn <- CrossTable(x = prey.testLabels, y = prey_pred, prop.chisq = FALSE)

print(prey_knn)

#determine if MARTEN/WEASEL prey groups are significantly different in C and N----

# Shapiro-Wilk normality test for C13 and N15 by prey group
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

with(dat.prey, shapiro.test(N[group_2 == "smallmamm.bird"]))

with(dat.prey, shapiro.test(C[group_2 == "smallmamm.bird"]))

with(dat.prey, shapiro.test(N[group_2 == "squirrel"]))

with(dat.prey, shapiro.test(C[group_2 == "squirrel"]))

with(dat.prey, shapiro.test(N[group_2 == "berry"]))

with(dat.prey, shapiro.test(C[group_2 == "berry"]))

#kruskal-wallis test of differences between > 2 groups in C13

dat.prey$group_2 <- as.factor(dat.prey$group_2)

kruskal.test(C ~ group_2, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

kwAllPairsNemenyiTest(data = dat.prey, x = C, g = group_2, dist="Tukey")

#               berry  smallmamm.bird
#smallmamm.bird 0.0290 -             
#squirrel       0.0011 0.0647      

#kruskal-wallis test of differences between > 2 groups in N15

dat.prey$group_2 <- as.factor(dat.prey$group_2)

kruskal.test(N ~ group_2, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

kwAllPairsNemenyiTest(data = dat.prey, x = N, g = group_2, dist="Tukey")

#               berry   smallmamm.bird
#smallmamm.bird 7.3e-06 -             
#squirrel       0.0022  0.6839     

#Run mixing models with MARTEN HAIR AND NAIL as consumer----

# load SI data (cleaned) from hair of consumers
#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

rm(list = ls())

dat.marten <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                           sql = "select * from file where `species` = 'marten'", sep = ",")

#make values in Tissue column either hair or nail

dat.marten <- dat.marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.marten <- dat.marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and marten

dat.means <- setNames(aggregate(list(dat.marten$C, dat.marten$N), 
                                by = list(dat.marten$ID, dat.marten$Tissue, dat.marten$spps_code, 
                                          dat.marten$species, dat.marten$sex), mean),
                      c("ID", "Tissue", "spps_code", "species", "sex", "C", "N"))

#subset consumer data by tissue type manually

marten_h <- subset(dat.means, Tissue == 'hair')

marten_n <- subset(dat.means, Tissue == 'nail')

#replace labels in respective dataframes

marten_h <- marten_h %>% 
  mutate(species = ifelse(as.character(species) == "marten", "marten_h", as.character(species)))

marten_n <- marten_n %>% 
  mutate(species = ifelse(as.character(species) == "marten", "marten_n", as.character(species)))

#fix spps_codes so that marten_h = 1, marten_n = 2

marten_h <- marten_h %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 2, 1, as.character(spps_code)))

#make new dataframe 

dat.marten.tis <- rbind(marten_h, marten_n) 

consumer.spps <- dat.marten.tis[c("spps_code", "C", "N")]

names(consumer.spps)[2] <- "d13C"

names(consumer.spps)[3] <- "d15N"

consumer.matrix <- data.matrix(consumer.spps, rownames.force = NA)

#Add prey data (source data)

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove fish, beaver, muskrat, grouse, waterfowl, salamanders, and toads from analysis

dat.prey <- dat.prey[-c(4:6, 22:74), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Squirrel", "squirrel", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

#make new data frame with Mean and SD for C and N by prey group
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "smallmamm.bird", 1, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "squirrel", 2, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "berry", 3, as.character(group_2)))

dat.prey = data.table(dat.prey)

dat.source = dat.prey[, list("Meand13C" = mean(C),
                             "SDd13C" = sd(C),
                             "Meand15N" = mean(N), 
                             "SDd15N" = sd(N)), 
                      by = c("group_2")]

#convert to matrix

source.matrix <- data.matrix(dat.source)

#make trophic enrichment fractor table with TEF values for 4 sources in this model
#1 = smallmamm.bird, 2 = squirrel, 3 = berry
#dC13 = +2 for terrestrial vertebrates and berries
#dN15 = +3 for terrestrial vertebrates and berries

dat.tef <- matrix(c(1, 2, 0.682, 3, 1.02, 
                    2, 2, 0.682, 3, 1.02, 
                    3, 2, 0.682, 3, 1.02), ncol=5,byrow=TRUE)

colnames(dat.tef) <- c("Source", "Meand13C","SDd13C","Meand15N", "SDd15N")

tef.df <- as.data.frame.matrix(dat.tef)

tef.matrix <- data.matrix(tef.df, rownames.force = NA)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

species.mod <- siarmcmcdirichletv4(consumer.matrix, source.matrix, tef.matrix, concdep = 0, 500000, 50000)

#view raw output data

out <- species.mod$output

#TABLE. summarize model output and write to csv----

melted <- melt(out)

summary.species <- ddply(melted, c("Var2"), summarise,
                         mean = mean(value), max = max(value), min = min(value), sd = sd(value),
                         CI = 1.96*(sd(value)/sqrt(length(value))), 
                         quant_5 = quantile(value, 0.05), 
                         quant_95 = quantile(value, 0.95), 
                         count = length(value))

write.csv(summary.species, "outputs/mixing model outputs/summary_marten_tissues_10nov2022.csv")

#ESM 3. plot density distributions for MARTEN HAIR AND NAIL----

martenh_mod <- melted[1:90000,1:3]

martenn_mod <- melted[150001:240000,1:3]

# Save an object to a file

save(species.mod, melted, martenh_mod, martenn_mod, file = "RData/tissue_marten_27july2021.RData")

# Restore the objects

load(file = "RData/tissue_marten_27july2021.RData")

# Rename the column and the values in the factor
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

levels(martenh_mod$Var2)[levels(martenh_mod$Var2)=="1G1"] <- "smallmamm.bird"
levels(martenh_mod$Var2)[levels(martenh_mod$Var2)=="2G1"] <- "squirrel"
levels(martenh_mod$Var2)[levels(martenh_mod$Var2)=="3G1"] <- "berry"
names(martenh_mod)[names(martenh_mod)=="Var2"]  <- "Diet_Groups"

levels(martenn_mod$Var2)[levels(martenn_mod$Var2)=="1G2"] <- "smallmamm.bird"
levels(martenn_mod$Var2)[levels(martenn_mod$Var2)=="2G2"] <- "squirrel"
levels(martenn_mod$Var2)[levels(martenn_mod$Var2)=="3G2"] <- "berry"
names(martenn_mod)[names(martenn_mod)=="Var2"]  <- "Diet_Groups"

martenh_diet <- ggplot(data=martenh_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=9, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=8),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=8),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=8))
martenh_diet <- martenh_diet + guides(fill=guide_legend(title="Diet Groups"))

martenn_diet <- ggplot(data=martenn_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=9, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=8),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=8),
        legend.position="none",
        legend.text = element_blank())

tiff(file = "figures/supplementary material/ESM_3_27july2021.tif", 
     units = "mm", width = 156, height = 234, res = 300)

ggarrange(martenh_diet, martenn_diet,
          labels = c("A", "B"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=10, face = "bold", color = "black"),
          ncol = 1, nrow = 2)

dev.off()

#Run mixing models with WEASEL as consumer----

# load SI data (cleaned) from hair of consumers
#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

#remove any objects from previous analysis

rm(list = ls())

dat.weasel <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                           sql = "select * from file where `species` = 'weasel'", sep = ",")

#correct species code = 1

dat.weasel <- dat.weasel %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 3, 1, as.character(spps_code)))

consumer.spps <- dat.weasel[c("spps_code", "C", "N")]

names(consumer.spps)[2] <- "d13C"

names(consumer.spps)[3] <- "d15N"

consumer.matrix <- data.matrix(consumer.spps, rownames.force = NA)

#Add prey data (source data)

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove fish, beaver, muskrat, grouse, waterfowl, salamanders, and toads from analysis

dat.prey <- dat.prey[-c(4:6, 22:74), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Squirrel", "squirrel", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

#make new data frame with Mean and SD for C and N by prey group
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "smallmamm.bird", 1, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "squirrel", 2, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "berry", 3, as.character(group_2)))

dat.prey = data.table(dat.prey)

dat.source = dat.prey[, list("Meand13C" = mean(C),
                             "SDd13C" = sd(C),
                             "Meand15N" = mean(N), 
                             "SDd15N" = sd(N)), 
                      by = c("group_2")]

#convert to matrix

source.matrix <- data.matrix(dat.source)

#make trophic enrichment fractor table with TEF values for 4 sources in this model
#1 = smallmamm.bird, 2 = squirrel, 3 = berry
#dC13 = +2 for terrestrial vertebrates and berries
#dN15 = +3 for terrestrial vertebrates and berries

dat.tef <- matrix(c(1, 2, 0.682, 3, 1.02, 
                    2, 2, 0.682, 3, 1.02, 
                    3, 2, 0.682, 3, 1.02), ncol=5,byrow=TRUE)

colnames(dat.tef) <- c("Source", "Meand13C","SDd13C","Meand15N", "SDd15N")

tef.df <- as.data.frame.matrix(dat.tef)

tef.matrix <- data.matrix(tef.df, rownames.force = NA)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

species.mod <- siarmcmcdirichletv4(consumer.matrix, source.matrix, tef.matrix, concdep = 0, 500000, 50000)

#view raw output data

out <- species.mod$output

#TABLE. summarize model output and write to csv----

melted <- melt(out)

summary.species <- ddply(melted, c("Var2"), summarise,
                         mean = mean(value), max = max(value), min = min(value), sd = sd(value),
                         CI = 1.96*(sd(value)/sqrt(length(value))), 
                         quant_5 = quantile(value, 0.05), 
                         quant_95 = quantile(value, 0.95), 
                         count = length(value))

write.csv(summary.species, "outputs/mixing model outputs/summary_weasel_10nov2022.csv")

#ESM 4. plot density distributions for WEASEL----

weasel_mod <- melted[1:90000,1:3]

# Save an object to a file

save(species.mod, melted, weasel_mod, file = "RData/weasel_27july2021.RData")

# Restore the objects

load(file = "RData/weasel_27july2021.RData")

# Rename the column and the values in the factor
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="1"] <- "smallmamm.bird"
levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="2"] <- "squirrel"
levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="3"] <- "berry"
names(weasel_mod)[names(weasel_mod)=="Var2"]  <- "Diet_Groups"

weasel_diet <- ggplot(data=weasel_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=9, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=8),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=8),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=8))
weasel_diet <- weasel_diet + guides(fill=guide_legend(title="Diet Groups"))

tiff(file = "figures/supplementary material/ESM_4_27july2021.tif", 
     units = "mm", width = 174, height = 131, res = 300)

weasel_diet

dev.off()

#K Nearest Neighbor Randomization Test to assign OTTER prey groups----

#from: https://www.datacamp.com/community/tutorials/machine-learning-in-r
#acccessed on April 6, 2020

#Add prey data and plot as scatter plot

rm(list = ls())

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove salmon, kokanee, beaver, salamanders, toads from analysis

dat.prey <- dat.prey[-c(1:23, 29:31, 44:46, 57:60, 67:70, 75:94), ]

#split up into otter prey groups:
#muskrats, fish.waterfowl, clam

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Clam", "clam", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Burbot", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Lake Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

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
#k = sqrt(N)/2, therefore k = 3

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

prey_knn <- CrossTable(x = prey.testLabels, y = prey_pred, prop.chisq = FALSE)

print(prey_knn)

#determine if OTTER prey groups are significantly different in C and N----
#muskrats, waterfowl, fish
# Shapiro-Wilk normality test for C13 and N15 by prey group

with(dat.prey, shapiro.test(N[group_2 == "fish.waterfowl"]))

with(dat.prey, shapiro.test(C[group_2 == "fish.waterfowl"]))

with(dat.prey, shapiro.test(N[group_2 == "muskrat"]))

with(dat.prey, shapiro.test(C[group_2 == "muskrat"]))

with(dat.prey, shapiro.test(N[group_2 == "clam"]))

with(dat.prey, shapiro.test(C[group_2 == "clam"]))

#kruskal-wallis test of differences between > 2 groups in C13

dat.prey$group_2 <- as.factor(dat.prey$group_2)

kruskal.test(C ~ group_2, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

kwAllPairsNemenyiTest(data = dat.prey, x = C, g = group_2, dist="Tukey")

#               clam    fish.waterfowl
#fish.waterfowl 0.00662 -             
#muskrat        0.00084 0.09201   

#kruskal-wallis test of differences between > 2 groups in N15

dat.prey$group_2 <- as.factor(dat.prey$group_2)

kruskal.test(N ~ group_2, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

kwAllPairsNemenyiTest(x = N, g = group_2, dist="Tukey")

#               clam   fish.waterfowl
#fish.waterfowl 0.0021 -             
#muskrat        0.9045 0.1237  

#Run mixing models with otter as grouping variable----

# load SI data (cleaned) from hair of consumers
#spps_code: mink = 1, marten = 2, weasel = 3, otter = 4

#remove any objects from previous analysis

rm(list = ls())

dat.hair <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

#subset consumer data by species manually

otter <- subset(dat.hair, species == 'otter')

#fix spps_codes so that otter = 1

otter <- otter %>% 
  mutate(spps_code = ifelse(as.character(spps_code) == 4, 1, as.character(spps_code)))

#make new dataframe with f_mink, m_mink, marten, and otter

consumer.spps <- otter[c("spps_code", "C", "N")]

names(consumer.spps)[2] <- "d13C"

names(consumer.spps)[3] <- "d15N"

consumer.matrix <- data.matrix(consumer.spps, rownames.force = NA)

#Add prey data (source data)

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove salmon, kokanee, beaver, salamanders, toads from analysis

dat.prey <- dat.prey[-c(1:23, 29:31, 44:46, 57:60, 67:70, 75:94), ]

#split up into otter prey groups:
#muskrats, fish.waterfowl, clam

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Clam", "mussel", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Burbot", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Lake Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

#make new data frame with Mean and SD for C and N by prey group
#convert group_2 to numeric: fish.waterfowl = 1, mussel = 2, muskrat = 3

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "fish.waterfowl", 1, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "mussel", 2, as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(group_2) == "muskrat", 3, as.character(group_2)))

dat.prey = data.table(dat.prey)

dat.source = dat.prey[, list("Meand13C" = mean(C),
                             "SDd13C" = sd(C),
                             "Meand15N" = mean(N), 
                             "SDd15N" = sd(N)), 
                      by = c("group_2")]

#convert to matrix

source.matrix <- data.matrix(dat.source)

#make trophic enrichment fractor table with TEF values for 4 sources in this model
# fish.waterfowl = 1, mussels = 2, muskrat = 3
#dC13 = +2 for terrestrial vertebrates and berries
#       +1 for fish, mussels
#dN15 = +3 for terrestrial vertebrates and berries
#       +2 for fish, mussels

dat.tef <- matrix(c(1, 1, 0.341, 2, 0.682, 
                    2, 1, 0.341, 2, 0.682,
                    3, 2, 0.682, 3, 1.023), ncol=5,byrow=TRUE)

colnames(dat.tef) <- c("Source", "Meand13C","SDd13C","Meand15N", "SDd15N")

tef.df <- as.data.frame.matrix(dat.tef)

tef.matrix <- data.matrix(tef.df, rownames.force = NA)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

species.mod <- siarmcmcdirichletv4(consumer.matrix, source.matrix, tef.matrix, concdep = 0, 500000, 50000)

#view raw output data

out <- species.mod$output

#TABLE. summarize model output and write to csv----

melted <- melt(out)

summary.species <- ddply(melted, c("Var2"), summarise,
                         mean = mean(value), max = max(value), min = min(value), sd = sd(value),
                         CI = 1.96*(sd(value)/sqrt(length(value))), 
                         quant_5 = quantile(value, 0.05), 
                         quant_95 = quantile(value, 0.95), 
                         count = length(value))

write.csv(summary.species, "outputs/mixing model outputs/summary_otter_10nov2022.csv")

#ESM 5. plot density distributions for otter----

#density plots for otter

otter_mod <- melted[1:90000,1:3]

# Save an object to a file

save(species.mod, melted, otter_mod, file = "RData/otter_27july2021.RData")

# Restore the objects

load(file = "RData/otter_27july2021.RData")

# Rename the column and the values in the factor
# fish.waterfowl = 1, fish2 = 2, clam = 3, muskrat = 4

levels(otter_mod$Var2)[levels(otter_mod$Var2)=="1"] <- "fish.waterfowl"
levels(otter_mod$Var2)[levels(otter_mod$Var2)=="2"] <- "mussel"
levels(otter_mod$Var2)[levels(otter_mod$Var2)=="3"] <- "muskrat"
names(otter_mod)[names(otter_mod)=="Var2"]  <- "Diet_Groups"

otter_diet <- ggplot(data=otter_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,30), expand = c(0, 0)) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=9, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=9),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=8),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=8),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=8))
otter_diet <- otter_diet + guides(fill=guide_legend(title="Diet Groups"))

tiff(file = "figures/supplementary material/ESM_5_27july2021.tif", 
     units = "mm", width = 174, height = 131, res = 300)

otter_diet

dev.off()

#Figure - Plot all density distributions in a single figure with 6 panels----

#remove any objects from previous analysis

rm(list = ls())

#Load in mixing model results for each species/tissue

load(file = "RData/tissue_mink_27july2021.RData")
load(file = "RData/tissue_marten_27july2021.RData")
load(file = "RData/weasel_27july2021.RData")
load(file = "RData/otter_27july2021.RData")

# Rename the column and the values in the factor
#convert group_2 to numeric: #berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4

levels(minkh_mod$Var2)[levels(minkh_mod$Var2)=="1G1"] <- "berry"
levels(minkh_mod$Var2)[levels(minkh_mod$Var2)=="2G1"] <- "fish.waterfowl"
levels(minkh_mod$Var2)[levels(minkh_mod$Var2)=="3G1"] <- "smallmamm.bird"
levels(minkh_mod$Var2)[levels(minkh_mod$Var2)=="4G1"] <- "muskrat"
names(minkh_mod)[names(minkh_mod)=="Var2"]  <- "Diet_Groups"

levels(minkn_mod$Var2)[levels(minkn_mod$Var2)=="1G2"] <- "berry"
levels(minkn_mod$Var2)[levels(minkn_mod$Var2)=="2G2"] <- "fish.waterfowl"
levels(minkn_mod$Var2)[levels(minkn_mod$Var2)=="3G2"] <- "smallmamm.bird"
levels(minkn_mod$Var2)[levels(minkn_mod$Var2)=="4G2"] <- "muskrat"
names(minkn_mod)[names(minkn_mod)=="Var2"]  <- "Diet_Groups"

# Rename the column and the values in the factor
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

levels(martenh_mod$Var2)[levels(martenh_mod$Var2)=="1G1"] <- "smallmamm.bird"
levels(martenh_mod$Var2)[levels(martenh_mod$Var2)=="2G1"] <- "squirrel"
levels(martenh_mod$Var2)[levels(martenh_mod$Var2)=="3G1"] <- "berry"
names(martenh_mod)[names(martenh_mod)=="Var2"]  <- "Diet_Groups"

levels(martenn_mod$Var2)[levels(martenn_mod$Var2)=="1G2"] <- "smallmamm.bird"
levels(martenn_mod$Var2)[levels(martenn_mod$Var2)=="2G2"] <- "squirrel"
levels(martenn_mod$Var2)[levels(martenn_mod$Var2)=="3G2"] <- "berry"
names(martenn_mod)[names(martenn_mod)=="Var2"]  <- "Diet_Groups"

# Rename the column and the values in the factor
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="1"] <- "smallmamm.bird"
levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="2"] <- "squirrel"
levels(weasel_mod$Var2)[levels(weasel_mod$Var2)=="3"] <- "berry"
names(weasel_mod)[names(weasel_mod)=="Var2"]  <- "Diet_Groups"

# Rename the column and the values in the factor
# fish.waterfowl = 1, mussel = 2, muskrat = 3

levels(otter_mod$Var2)[levels(otter_mod$Var2)=="1"] <- "fish.waterfowl"
levels(otter_mod$Var2)[levels(otter_mod$Var2)=="2"] <- "mussel"
levels(otter_mod$Var2)[levels(otter_mod$Var2)=="3"] <- "muskrat"
names(otter_mod)[names(otter_mod)=="Var2"]  <- "Diet_Groups"

#Plot the results

#Diet Group Fill Colors:
#1. Berry = red
#2. Fish.waterfowl = blue
#3. smallmamm.bird = green
#4. muskrat = orange
#5. squirrel = purple
#6. mussel = pink

#Mink Hair (A)
#berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4

minkh_diet <- ggplot(data=minkh_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0, 0)) +
  scale_fill_manual(values = c("green","orange", "blue", "red" )) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=10, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=10),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=9),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=9),
        legend.position = c(0.9, 0.9),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=9))
minkh_diet <- minkh_diet + guides(fill=guide_legend(title="Diet Groups"))

#Marten hair (B)
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

martenh_diet <- ggplot(data=martenh_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0, 0)) +
  scale_fill_manual(values = c("green","purple", "red")) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=10, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=10),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=9),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=9),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size= 9))
martenh_diet <- martenh_diet + guides(fill=guide_legend(title="Diet Groups"))

#Weasel (C)
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

weasel_diet <- ggplot(data=weasel_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0, 0)) +
  scale_fill_manual(values = c("green","purple", "red")) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=10, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=10),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=9),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=9),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=9))
weasel_diet <- weasel_diet + guides(fill=guide_legend(title="Diet Groups"))

#Mink Nail (D)
#berry = 1, fish.waterfowl = 2, smallmamm.bird = 3, muskrat = 4

minkn_diet <- ggplot(data=minkn_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 50), expand = c(0, 0)) +
  scale_fill_manual(values = c("green","orange", "blue", "red")) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=10, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=10),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=9),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=9),
        legend.position = c(0.9, 0.9),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=9))
minkn_diet <- minkn_diet + guides(fill=guide_legend(title="Diet Groups"))

#Marten nail (E)
#1 = smallmamm.bird, 2 = squirrel, 3 = berry

martenn_diet <- ggplot(data=martenn_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0, 0)) +
  scale_fill_manual(values = c("green","purple", "red")) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=10, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=10),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=9),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=9),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=9))
martenn_diet <- martenn_diet + guides(fill=guide_legend(title="Diet Groups"))

#Otter (F)
# fish.waterfowl = 1, mussel = 2, muskrat = 3

otter_diet <- ggplot(data=otter_mod, aes(x= value, group = Diet_Groups, fill = Diet_Groups)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,50), expand = c(0, 0)) +
  scale_fill_manual(values = c("orange","blue", "pink")) +
  labs(x = "Dietary Composition", y = "Density") +
  theme(axis.title.x = element_text(family="Times New Roman", size=10, face = "plain"),
        axis.title.y = element_text(family="Times New Roman", size=10),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(family = "Times New Roman", colour = "black", size=9),
        axis.text.y = element_text(family = "Times New Roman", colour = "black", size=9),
        legend.position = c(0.9, 0.92),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = 9, face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size=9))
otter_diet <- otter_diet + guides(fill=guide_legend(title="Diet Groups"))

tiff(file = "figures/supplementary material/Figs_5_thru_8_20230726.tif", 
     units = "mm", width = 356, height = 234, res = 300)

ggarrange(minkh_diet, martenh_diet, weasel_diet, minkn_diet, martenn_diet, otter_diet,
          labels = c("A. Mink - hair", "C. Marten - hair", "E. Weasel - hair", "B. Mink - nail",
                     "D. Marten - nail", "F. Otter - hair"),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=10, face = "bold", color = "black"),
          ncol = 3, nrow = 2)

dev.off()

#Figure - Plot Fig 2 thru 4 as single panel----

#Figure 2

rm(list = ls())

#load SI data (cleaned) from consumers

dat.mustel <- read.csv("data/jprf_consumer_data_7jun2021.csv", header = TRUE, sep = ",")

#make values in Tissue column either hair or nail

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mustel <- dat.mustel %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and marten

dat.means <- setNames(aggregate(list(dat.mustel$C, dat.mustel$N), 
                                by = list(dat.mustel$ID, dat.mustel$Tissue, dat.mustel$spps_code, 
                                          dat.mustel$species), mean),
                      c("ID", "Tissue", "spps_code", "species", "C", "N"))

#subset consumer data by species manually

mink_h <- subset(dat.means, species == 'mink' & Tissue == 'hair')

mink_n <- subset(dat.means, species == 'mink' & Tissue == 'nail')

marten <- subset(dat.means, species == 'marten'& Tissue == 'hair')

weasel <- subset(dat.means, species == 'weasel')

otter <- subset(dat.means, species == 'otter')

#make new dataframe with all mustelid hair data

dat.mustel.hair <- rbind(mink_h, marten, weasel, otter) 

#Figure 1A. Join means by tissues for mustelids to the original data frame and pipe to ggplot

consumer_means <- left_join(dat.mustel.hair, 
                            dat.mustel.hair %>%
                              group_by(species) %>%
                              summarise_at(vars(C, N), funs(mean = mean))
) %>% 
  ggplot(aes(colour=species, shape = species)) +
  # Plot raw data
  geom_point(aes(x=C, y=N), size = 3, alpha = 0.8) +
  # Plot mean by species
  geom_point(aes(x=C_mean, y=N_mean), size=4) +
  # Add vertical errorbars 
  scale_shape_manual(values = c(15, 17, 18, 19)) +
  scale_color_manual(values = c("#D55E00", "#0072B2", "#009E73", "#E69F00")) + 
  stat_summary(fun.data = mean_se, fun.args = list(mult=1.96), 
               aes(x = C_mean, y=N),
               geom="errorbar", width=0.1) +
  # Add horizontal errorbars
  stat_summaryh(fun.data=mean_se_h, fun.args=list(mult=1.96), 
                aes(x=C, y=N_mean),
                geom="errorbarh", width=0.1) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  #labs(colour = "Mustelid Groups", shape = "Mustelid Groups") +
  #labs(title = "Isotopic Signatures", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=15, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=14, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        panel.background = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "15",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.key = element_blank(),
        legend.position = c(0.90, 0.95),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_continuous(breaks = seq(-32, -18, 1), limits = c(-32,-18), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0, 15, 1), limits = c(-0,15), expand = c(0, 0))

consumer_means <- consumer_means + labs(fill = "Mustelid Species", 
                                        colour = "Mustelid Species", shape = "Mustelid Species")

#Figure 3 - ellipses with prey values

rm(dat.means, dat.mustel, dat.mustel.hair, marten, mink_h, mink_n, otter, weasel)

dat.hair <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `Tissue` = 'hair'", sep = ",")

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

SE1
SE2
SE3
SE4

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

mink_ <- rep("mink", length(SE1$xSEAc))
marten_ <- rep("marten", length(SE2$xSEAc))
weasel_ <- rep("weasel", length(SE3$xSEAc))
otter_ <- rep("otter", length(SE4$xSEAc))

#create new data frame with names and ellipse outputs

species <- c(mink_, marten_, weasel_, otter_)
x <- c(SE1$xSEAc,SE2$xSEAc, SE3$xSEAc, SE4$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc, SE3$ySEAc, SE4$ySEAc)

df_SE <- data.frame(x, y, species)

plot(df_SE$x, df_SE$y)

#plot SE1$xSEAc vs SE1$ySEAc; the ellipse(s) form as a series of points
#This gives SEAc (sample size corrected area)
#hull total area (TA)

find_hull <- function(df) df[chull(df$C, df$N), ]

hulls <- ddply(dat.hair, "species", find_hull)

#Add prey data and plot 

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, salamanders, toads, lake trout, salmon, kokanee, squirrel, grouse from analysis

dat.prey <- dat.prey[-c(4:6, 22:23, 26:31, 44:49, 57:60, 67:70), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Squirrel", "squirrel", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Clam", "mussel", as.character(group_2)))

#correct prey by TEF values as in mixing models
#dC13 = +2 for terrestrial vertebrates and berries
#       +1 for fish, clams
#dN15 = +3 for terrestrial vertebrates and berries
#       +2 for fish, clams

dat.prey$C <- dat.prey$C + ifelse(dat.prey$group_2 == "fish.waterfowl", 1, 
                                  ifelse(dat.prey$group_2 == "mussel", 1, 2))

dat.prey$N <- dat.prey$N + ifelse(dat.prey$group_2 == "fish.waterfowl", 2, 
                                  ifelse(dat.prey$group_2 == "mussel", 2, 3))

#Combine the mustelid ellipses with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(group_2) %>%
                             dplyr::summarize(mean_C = mean(C),
                                              se_C = sd(C)/sqrt(n()),
                                              mean_N = mean(N),
                                              se_N=sd(N)/sqrt(n()),
                                              sample_size = n()))

mustelid_w_prey <- ggplot(NULL) +
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = group_2), size = 4) +
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_path(data = df_SE, aes(x = x, y = y, color = species), linetype = 1, size = 1) +
  
  scale_shape_manual(values = c(1,2,0,5,6,7)) +
  scale_color_manual(values = c("#D55E00", "#0072B2", "#009E73", "#E69F00")) + 
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  #labs(title = "Mustelid Hair with Prey", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=15, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=14, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.90, 0.89), legend.direction = "vertical",
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "15",face = "bold"),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2)) +
  scale_x_continuous(breaks = seq(-32, -18, 1), limits = c(-32,-18), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0, 15, 1), limits = c(-0,15), expand = c(0, 0))

mustelid_w_prey <- mustelid_w_prey + labs(shape = "Prey", color = "Mustelid Species")

#Figure 4 - marten and mink hair and nail with prey means

rm(dat.hair, dat.prey, df_SE, hulls, marten, mink, otter, prey_mean, SE1, SE2, SE3, SE4, weasel,
   marten_, mink_, otter_, species, weasel_, x, y, find_hull)

dat.mima <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `species` = 'mink' 
                         or `species` = 'marten'", sep = ",")

#make values in Tissue column either hair or nail

dat.mima <- dat.mima %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "nail", as.character(Tissue)))

dat.mima <- dat.mima %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "nail", as.character(Tissue)))

#Take mean of bottom and top nail segments for mink and marten

dat.mima <- setNames(aggregate(list(dat.mima$C, dat.mima$N), 
                               by = list(dat.mima$ID, dat.mima$Tissue, dat.mima$spps_code, 
                                         dat.mima$species, dat.mima$sex), mean),
                     c("ID", "Tissue", "spps_code", "species", "sex", "C", "N"))

#subset consumer data by species (mink and marten) manually

mink <- subset(dat.mima, species == 'mink')

marten <- subset(dat.mima, species == 'marten')

#subset consumer data by tissue type manually

mink.hair <- subset(mink, Tissue == 'hair')

mink.nail <- subset(mink, Tissue == 'nail')

marten.hair <- subset(marten, Tissue == 'hair')

marten.nail <- subset(marten, Tissue == 'nail')

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(mink.hair$C, mink.hair$N)
SE2 <- standard.ellipse(mink.nail$C, mink.nail$N)
SE3 <- standard.ellipse(marten.hair$C, marten.hair$N)
SE4 <- standard.ellipse(marten.nail$C, marten.nail$N)

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

mink.hair_ <- rep("mink.hair", length(SE1$xSEAc))
mink.nail_ <- rep("mink.nail", length(SE2$xSEAc))
marten.hair_ <- rep("marten.hair", length(SE3$xSEAc))
marten.nail_ <- rep("marten.nail", length(SE4$xSEAc))

#create new data frame with names and ellipse outputs

tissues <- c(mink.hair_, mink.nail_, marten.hair_, marten.nail_)
x <- c(SE1$xSEAc,SE2$xSEAc, SE3$xSEAc, SE4$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc,SE3$ySEAc)

df_SE <- data.frame(x, y, tissues)

plot(df_SE$x, df_SE$y)

#Add prey data and plot 

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, clam, salamanders, toads, lake trout, salmon, kokanee, grouse from analysis

dat.prey <- dat.prey[-c(4:6, 22:23, 26:31, 35:38, 44:49, 57:60, 67:70), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Squirrel", "squirrel", as.character(group_2)))

#correct prey by TEF values as in mixing models
#dC13 = +2 for terrestrial vertebrates and berries
#       +1 for fish
#dN15 = +3 for terrestrial vertebrates and berries
#       +2 for fish

dat.prey$C <- dat.prey$C + ifelse(dat.prey$group_2 == "fish.waterfowl", 1, 2)

dat.prey$N <- dat.prey$N + ifelse(dat.prey$group_2 == "fish.waterfowl", 2, 3)

#Combine the mink ellipses and hulls with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(group_2) %>%
                             dplyr::summarize(mean_C = mean(C),
                                              se_C = sd(C)/sqrt(n()),
                                              mean_N = mean(N),
                                              se_N=sd(N)/sqrt(n()),
                                              sample_size = n()))

mima_tissues_w_prey <- ggplot(NULL) +
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = group_2), size = 4) +
  
  
  geom_path(data = df_SE, aes(x = x, y = y, color = tissues), linetype = 1, size = 1) +
  
  scale_shape_manual(values = c(1,0,5,6,7)) +
  scale_color_manual(values = c("#D55E00", "#CC79A7", "#0072B2", "#56B4E9")) +
  
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  #labs(title = "Mink and Marten Tissues with Prey", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=15, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=14, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.90, 0.89), legend.direction = "vertical",
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "15",face = "bold"),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.spacing = unit(0, 'cm'),
        legend.key.size = unit(0.3, "cm")) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2)) +
  scale_x_continuous(breaks = seq(-32, -18, 1), limits = c(-32,-18), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0, 15, 1), limits = c(-0,15), expand = c(0, 0))

mima_tissues_w_prey <- mima_tissues_w_prey + labs(shape = "Prey", color = "Mink and Marten Tissues")

tiff(file = "figures/supplementary material/Figs_2_thru_4_20230728.tif", 
     units = "mm", width = 1056, height = 334, res = 300)

ggarrange(consumer_means, mustelid_w_prey, mima_tissues_w_prey,
          labels = c("A.", "B.", "C."),
          hjust = -1,
          vjust = 2,
          font.label = list(family="Times New Roman", size=15, face = "bold", color = "black"),
          ncol = 3, nrow = 1)

dev.off()

#FIGURE Plot OTTER ellipses and convex hulls with prey (mean +/- SE)----

rm(list = ls())

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

#subset consumer data by species manually

otter <- subset(dat.hair, species == 'otter')

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(otter$C, otter$N)

SE1

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

otter_ <- rep("otter", length(SE1$xSEAc))

#create new data frame with names and ellipse outputs

species <- c(otter_)
x <- c(SE1$xSEAc)
y <- c(SE1$ySEAc)

df_SE <- data.frame(x, y, species)

plot(df_SE$x, df_SE$y)

#plot SE1$xSEAc vs SE1$ySEAc; the ellipse(s) form as a series of points
#This gives SEAc (sample size corrected area)
#hull total area (TA)

find_hull <- function(df) df[chull(df$C, df$N), ]

hulls <- ddply(otter, "species", find_hull)

#Add prey data and plot 

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove salmon, kokanee, beaver, salamanders, toads from analysis

dat.prey <- dat.prey[-c(1:23, 29:31, 44:46, 57:60, 67:70, 75:94), ]

#split up into otter prey groups:
#muskrats, kokanee, salmon, aquatic.prey

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish1", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish1", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish1", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Clam", "fish1", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish1", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish1", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish1", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Burbot", "fish2", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Lake Trout", "fish2", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "waterfowl", as.character(group_2)))

#correct prey by TEF values as in mixing models
#dC13 = +2 for terrestrial vertebrates and berries
#       +1 for fish, clams
#dN15 = +3 for terrestrial vertebrates and berries
#       +2 for fish, clams

dat.prey$C <- dat.prey$C + ifelse(dat.prey$group_2 == "fish1", 1, 
                                  ifelse(dat.prey$group_2 == "fish2", 1, 2))

dat.prey$N <- dat.prey$N + ifelse(dat.prey$group_2 == "fish1", 2, 
                                  ifelse(dat.prey$group_2 == "fish2", 2, 3))

#Combine the mustelid ellipses and hulls with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(group_2) %>%
                             dplyr::summarize(mean_C = mean(C),
                                              se_C = sd(C)/sqrt(n()),
                                              mean_N = mean(N),
                                              se_N=sd(N)/sqrt(n()),
                                              sample_size = n()))

otter_w_prey <- ggplot(NULL) +
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = group_2), size = 4) +
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_polygon(data = hulls, aes(x = C, y = N, color = species), alpha = 0.1, linetype = 3) +
  geom_path(data = df_SE, aes(x = x, y = y, color = species), linetype = 1, size = 1) +
  
  scale_shape_manual(values = c(1,2,0,5)) +
  scale_color_manual(values = c("blue", "red", "green", "purple", "orange")) + 
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Otter Hair with Prey", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(1.1, .85), legend.direction = "vertical",
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "16",face = "bold"),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2)) +
  scale_x_continuous(breaks = seq(-33, -18, 1)) +
  scale_y_continuous(breaks = seq(-0, 13, 1))

otter_w_prey <- otter_w_prey + labs(shape = "Prey", color = "River Otter")

tiff(file = "figures/hulls with prey/otter_prey_corr_16jun2020.tif", 
     units = "cm", width = 25, height = 20, res = 300)

otter_w_prey

dev.off()

#Plot ALL prey data as scatter with LABELS for REFERENCE----
#Add prey data and plot as scatter plot
#remove any objects from previous analysis

rm(list = ls())

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, salamanders, toads, and lake trout from analysis

dat.prey <- dat.prey[-c(22:23, 29:31, 44:46, 57:60, 67:70), ]

#plot prey by species

prey_species <- ggplot(dat.prey, aes(x = C , y = N)) +
  geom_point(aes(shape = species, color = species)) +
  geom_label_repel(aes(label = species),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  scale_x_continuous(breaks = seq(-34, -18, 1)) +
  scale_y_continuous(breaks = seq(-5, 13, 1)) +
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
        legend.position="none",
        legend.text = element_blank())

#Export plot

tiff(file = "figures/prey_species_17jun2020.tif", 
     units = "cm", width = 30, height = 17, res = 300)

prey_species

dev.off()

#Figure. Plot Mink and Marten Hair and Nail ellipses BY SEX with prey (Mean +/- SE)----

rm(list = ls())

dat.mima <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `species` = 'mink' 
                         or `species` = 'marten'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.mima$C, plot = FALSE)$out

dat.mima[which(dat.mima$C %in% outlier),]

dat.mima <- dat.mima[-which(dat.mima$C %in% outlier),]

outlier <- boxplot(dat.mima$N, plot = FALSE)$out

dat.mima[which(dat.mima$N %in% outlier),]

dat.mima <- dat.mima[-which(dat.mima$N %in% outlier),]

boxplot(dat.mima$C)

boxplot(dat.mima$N)

#subset consumer data by species (mink and marten) manually

mink <- subset(dat.mima, species == 'mink')

marten <- subset(dat.mima, species == 'marten')

#make values in Tissue column either hair or nail

mink <- mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "mink.nail", as.character(Tissue)))

mink <- mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "mink.nail", as.character(Tissue)))

mink <- mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "hair", "mink.hair", as.character(Tissue)))

marten <- marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "marten.nail", as.character(Tissue)))

marten <- marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "marten.nail", as.character(Tissue)))

marten <- marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "hair", "marten.hair", as.character(Tissue)))

#subset consumer data by tissue type manually

mink.hair <- subset(mink, Tissue == 'mink.hair')

mink.nail <- subset(mink, Tissue == 'mink.nail')

marten.hair <- subset(marten, Tissue == 'marten.hair')

marten.nail <- subset(marten, Tissue == 'marten.nail')

f.mink.hair <- subset(mink.hair, sex == 'F')

m.mink.hair <- subset(mink.hair, sex == 'M')

f.mink.nail <- subset(mink.nail, sex == 'F')

m.mink.nail <- subset(mink.nail, sex == 'M')

f.marten.hair <- subset(marten.hair, sex == 'F')

m.marten.hair <- subset(marten.hair, sex == 'M')

f.marten.nail <- subset(marten.nail, sex == 'F')

m.marten.nail <- subset(marten.nail, sex == 'M')


#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(mink.hair$C, mink.hair$N)
SE2 <- standard.ellipse(mink.nail$C, mink.nail$N)
SE3 <- standard.ellipse(marten.hair$C, marten.hair$N)
SE4 <- standard.ellipse(marten.nail$C, marten.nail$N)

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

mink.hair_ <- rep("mink.hair", length(SE1$xSEAc))
mink.nail_ <- rep("mink.nail", length(SE2$xSEAc))
marten.hair_ <- rep("marten.hair", length(SE3$xSEAc))
marten.nail_ <- rep("marten.nail", length(SE4$xSEAc))

#create new data frame with names and ellipse outputs

tissues <- c(mink.hair_, mink.nail_, marten.hair_, marten.nail_)
x <- c(SE1$xSEAc,SE2$xSEAc, SE3$xSEAc, SE4$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc,SE3$ySEAc)

df_SE <- data.frame(x, y, tissues)

plot(df_SE$x, df_SE$y)

#Add prey data and plot 

dat.prey <- read.csv("data/jprf_prey_data_17apr2020.csv", header = TRUE, sep = ",")

#remove burbot, beaver, clam, salamanders, toads, lake trout, salmon, kokanee, grouse from analysis

dat.prey <- dat.prey[-c(4:6, 22:23, 26:31, 35:38, 44:49, 57:60, 67:70), ]

#split terrestrial vertebrates into aquatic vs terrestrial vertebrates

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Goose", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Grebe", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Teal", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Northern Pikeminnow", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Rainbow Trout", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Chub", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sculpin", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Sucker", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Whitefish", "fish.waterfowl", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Vole", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Mouse", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Hare", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Shrew", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "SongBird", "smallmamm.bird", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Muskrat", "muskrat", as.character(group_2)))

dat.prey <- dat.prey %>% 
  mutate(group_2 = ifelse(as.character(species) == "Squirrel", "squirrel", as.character(group_2)))

#correct prey by TEF values as in mixing models
#dC13 = +2 for terrestrial vertebrates and berries
#       +1 for fish
#dN15 = +3 for terrestrial vertebrates and berries
#       +2 for fish

dat.prey$C <- dat.prey$C + ifelse(dat.prey$group_2 == "fish.waterfowl", 1, 2)

dat.prey$N <- dat.prey$N + ifelse(dat.prey$group_2 == "fish.waterfowl", 2, 3)

#Combine the mink ellipses and hulls with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(group_2) %>%
                             dplyr::summarize(mean_C = mean(C),
                                              se_C = sd(C)/sqrt(n()),
                                              mean_N = mean(N),
                                              se_N=sd(N)/sqrt(n()),
                                              sample_size = n()))

mima_tissues_w_prey <- ggplot(NULL) +
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = group_2), size = 4) +
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_path(data = df_SE, aes(x = x, y = y, color = tissues), linetype = 1, size = 1) +
  
  scale_shape_manual(values = c(1,0,5,6,7)) +
  scale_color_manual(values = c("blue", "red", "green", "yellow")) + 
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Mink and Marten Tissues with Prey", subtitle = "mean ± 95% CI") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.9, .3), legend.direction = "vertical",
        legend.title = element_text(family = "Times New Roman", colour = "black", size = "16",face = "bold"),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2)) +
  scale_x_continuous(breaks = seq(-33, -18, 1)) +
  scale_y_continuous(breaks = seq(-0, 13, 1))

mima_tissues_w_prey <- mima_tissues_w_prey + labs(shape = "Prey", color = "Mink and Marten Tissues")

tiff(file = "figures/figure3_08jun2021.tif", 
     units = "cm", width = 25, height = 20, res = 300)

mima_tissues_w_prey

dev.off()

#TABLE. calculate % niche overlap between consumers incl. marten nail and hair by SEX----

rm(list = ls())

dat.mima <- read.csv.sql("data/jprf_consumer_data_7jun2021.csv", header = TRUE, 
                         sql = "select * from file where `species` = 'mink' 
                         or `species` = 'marten'", sep = ",")

#Find and remove outlier from data set (outlier marten C value = MT138)

outlier <- boxplot(dat.mima$C, plot = FALSE)$out

dat.mima[which(dat.mima$C %in% outlier),]

dat.mima <- dat.mima[-which(dat.mima$C %in% outlier),]

outlier <- boxplot(dat.mima$N, plot = FALSE)$out

dat.mima[which(dat.mima$N %in% outlier),]

dat.mima <- dat.mima[-which(dat.mima$N %in% outlier),]

boxplot(dat.mima$C)

boxplot(dat.mima$N)

#subset consumer data by species (mink and marten) manually

mink <- subset(dat.mima, species == 'mink')

marten <- subset(dat.mima, species == 'marten')

#make values in Tissue column either hair or nail

mink <- mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "mink.nail", as.character(Tissue)))

mink <- mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "mink.nail", as.character(Tissue)))

mink <- mink %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "hair", "mink.hair", as.character(Tissue)))

marten <- marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailb", "marten.nail", as.character(Tissue)))

marten <- marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "nailt", "marten.nail", as.character(Tissue)))

marten <- marten %>% 
  mutate(Tissue = ifelse(as.character(Tissue) == "hair", "marten.hair", as.character(Tissue)))

#subset consumer data by tissue type manually

mink.hair <- subset(mink, Tissue == 'mink.hair')

mink.nail <- subset(mink, Tissue == 'mink.nail')

marten.hair <- subset(marten, Tissue == 'marten.hair')

marten.nail <- subset(marten, Tissue == 'marten.nail')

f.mink.hair <- subset(mink.hair, sex == 'F')

m.mink.hair <- subset(mink.hair, sex == 'M')

f.mink.nail <- subset(mink.nail, sex == 'F')

m.mink.nail <- subset(mink.nail, sex == 'M')

f.marten.hair <- subset(marten.hair, sex == 'F')

m.marten.hair <- subset(marten.hair, sex == 'M')

f.marten.nail <- subset(marten.nail, sex == 'F')

m.marten.nail <- subset(marten.nail, sex == 'M')

#make values in Tissue column either hair or nail

f.mink.hair <- data.frame(append(f.mink.hair, c(tissue_sex='f.mink.hair')))

m.mink.hair <- data.frame(append(m.mink.hair, c(tissue_sex='m.mink.hair')))

f.mink.nail <- data.frame(append(f.mink.nail, c(tissue_sex='f.mink.nail')))

m.mink.nail <- data.frame(append(m.mink.nail, c(tissue_sex='m.mink.nail')))

f.marten.hair <- data.frame(append(f.marten.hair, c(tissue_sex='f.marten.hair')))

m.marten.hair <- data.frame(append(m.marten.hair, c(tissue_sex='m.marten.hair')))

f.marten.nail <- data.frame(append(f.marten.nail, c(tissue_sex='f.marten.nail')))

m.marten.nail <- data.frame(append(m.marten.nail, c(tissue_sex='m.marten.nail')))

mink.hair <- data.frame(append(mink.hair, c(tissue_sex='mink.hair')))

mink.nail <- data.frame(append(mink.nail, c(tissue_sex='mink.nail')))

marten.hair <- data.frame(append(marten.hair, c(tissue_sex='marten.hair')))

marten.nail <- data.frame(append(marten.nail, c(tissue_sex='marten.nail')))

#make new dataframe with all mustelid groups for comparison

dat.mustel.all <- rbind(f.mink.hair, m.mink.hair, f.mink.nail, m.mink.nail, 
                        f.marten.hair, m.marten.hair, f.marten.nail, m.marten.nail, 
                        mink.hair, mink.nail, marten.hair, marten.nail)

aggregate(dat.mustel.all[8:9], dat.mustel.all[10], mean) # isotope means per sex/tissue

# format data for plotting function

dat.mustel.all$tissue_sex <- as.factor(dat.mustel.all$tissue_sex)

# niche overlap plots for 95% niche region sizes
# generate parameter draws from the "default" posteriors of each species

nsamples <- 500

mustel.par <- tapply(1:nrow(dat.mustel.all), dat.mustel.all$tissue_sex,
                       function(ii) niw.post(nsamples = nsamples, X = dat.mustel.all[ii,8:9]))

# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.

over <- overlap(mustel.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, .99))

# posterior expectations of overlap metrics

over.mean <- apply(over*100, c(1:2, 4), mean)

round(over.mean)

write.csv(over.mean, "outputs/over.mean_tissuesex_08jun2021.csv")
