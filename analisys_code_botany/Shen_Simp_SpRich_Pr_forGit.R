library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggalt)
library(ggthemr)

setwd("/Users/admin/Yandex.Disk.localized/AnalisysSoromotin/")
df <- read.csv("GeoBotany_fic.sp_newName_newLocacion.csv", row.names = 1)
# loading a dataframe with initial data

write.csv(simp, 'Simpson.csv')
write.csv(shan, 'Shannon.csv') # writing dataframes to text files

## Biodiversity indices

# Simpson
simpson <- diversity(df, index = "simpson")
simp <- as.data.frame(simpson)
simp$site <- row.names(simp)
simp$site <- gsub("DuneTop(\\d|\\d{2})", "DuneTop", simp$site)
simp$site <- gsub("LeeSlo(\\d|\\d{2})", "LeeSlo", simp$site)
simp$site <- gsub("InterLow(\\d|\\d{2})", "InterLow", simp$site)
simp$site <- gsub("WindSlo(\\d|\\d{2})", "WindSlo", simp$site)
simp$site <- gsub("YoungF(\\d|\\d{2})", "YoungF", simp$site)
simp$site <- gsub("MatureF(\\d|\\d{2})", "MatureF", simp$site)
simp$site <- gsub("ClimaxF(\\d|\\d{2})", "ClimaxF", simp$site)

simp$site <- as.factor(simp$site)
simp <- mutate(simp, site = factor(
  site, levels = c("DuneTop", "LeeSlo","InterLow", "WindSlo", "YoungF", 
                   "MatureF", "ClimaxF")))
ggthemr('solarized')
ggthemr_reset()
simp_plt <- 
  ggplot(data = simp, aes(x = site, y = simpson)) +
  #geom_violin(scale = "count", aes(fill = site), 
  #           size = 0.3, width = 1.8) +
  geom_boxplot(width = 0.4, fill = NA, size = 0.3, 
               weight = 1) +
  labs(fill="site", shape="site", colour="site") +
  labs(x = "site", y = "simpson") +
  ggtitle("Simpson index") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

plot(simp_plt)
#stat_summary(fun = "median", geom = "point")
ggsave(filename = "simp_plt_4.png", plot = simp_plt, dpi = 600, units = "in", 
       width = 7, height = 5)

# Shannon
shannon <- diversity(df, index = "shannon")
shan <- as.data.frame(shannon)
shan$site <- row.names(shan)
shan$site <- gsub("DuneTop(\\d|\\d{2})", "DuneTop", shan$site)
shan$site <- gsub("LeeSlo(\\d|\\d{2})", "LeeSlo", shan$site)
shan$site <- gsub("InterLow(\\d|\\d{2})", "InterLow", shan$site)
shan$site <- gsub("WindSlo(\\d|\\d{2})", "WindSlo", shan$site)
shan$site <- gsub("YoungF(\\d|\\d{2})", "YoungF", shan$site)
shan$site <- gsub("MatureF(\\d|\\d{2})", "MatureF", shan$site)
shan$site <- gsub("ClimaxF(\\d|\\d{2})", "ClimaxF", shan$site)

shan$site <- as.factor(shan$site)
shan <- mutate(shan, site = factor(
  site, levels = c("DuneTop", "LeeSlo","InterLow", "WindSlo", "YoungF", 
                   "MatureF", "ClimaxF")))
ggthemr('solarized')
ggthemr_reset()
shan_plt <- 
  ggplot(data = shan, aes(x = site, y = shannon)) +
  #geom_violin(scale = "count", aes(fill = site), 
  #           size = 0.3, width = 1.8) +
  geom_boxplot(width = 0.4, fill = NA, size = 0.3, 
               weight = 1) +
  labs(fill="site", shape="site", colour="site") +
  labs(x = "site", y = "shannon") +
  ggtitle("Shannon index") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

plot(shan_plt)
#stat_summary(fun = "median", geom = "point")
ggsave(filename = "shan_plt_4.png", plot = shan_plt, dpi = 600, units = "in", 
       width = 7, height = 5)


#Species richness

richness <- specnumber(df)
rich <- as.data.frame(richness)
rich$site <- row.names(rich)
rich$site <- gsub("DuneTop(\\d|\\d{2})", "DuneTop", rich$site)
rich$site <- gsub("LeeSlo(\\d|\\d{2})", "LeeSlo", rich$site)
rich$site <- gsub("InterLow(\\d|\\d{2})", "InterLow", rich$site)
rich$site <- gsub("WindSlo(\\d|\\d{2})", "WindSlo", rich$site)
rich$site <- gsub("YoungF(\\d|\\d{2})", "YoungF", rich$site)
rich$site <- gsub("MatureF(\\d|\\d{2})", "MatureF", rich$site)
rich$site <- gsub("ClimaxF(\\d|\\d{2})", "ClimaxF", rich$site)

rich$site <- as.factor(rich$site)
rich <- mutate(rich, site = factor(
  site, levels = c("DuneTop", "LeeSlo","InterLow", "WindSlo", "YoungF", 
                   "MatureF", "ClimaxF")))
ggthemr('solarized')
ggthemr_reset()
rich_plt <- 
  ggplot(data = rich, aes(x = site, y = richness)) +
  #geom_violin(scale = "count", aes(fill = site), 
  #           size = 0.3, width = 1.8) +
  geom_boxplot(width = 0.4, fill = NA, size = 0.3, 
               weight = 1) +
  labs(fill="site", shape="site", colour="site") +
  labs(x = "site", y = "Species richness") +
  ggtitle("Species richness") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

plot(rich_plt)
#stat_summary(fun = "median", geom = "point")
ggsave(filename = "rich_plt_4.png", plot = rich_plt, dpi = 600, units = "in", 
       width = 7, height = 5)

##Projective coverage

#Grass
dfG <- read.csv("Coverage - Grass.csv", row.names = 1)

projection <- dfG[,1]
pr <- dfG

pr$site <- row.names(pr)
pr$site <- gsub("DuneTop(\\d|\\d{2})", "DuneTop", pr$site)
pr$site <- gsub("LeeSlo(\\d|\\d{2})", "LeeSlo", pr$site)
pr$site <- gsub("InterLow(\\d|\\d{2})", "InterLow", pr$site)
pr$site <- gsub("WindSlo(\\d|\\d{2})", "WindSlo", pr$site)
pr$site <- gsub("YoungF(\\d|\\d{2})", "YoungF", pr$site)
pr$site <- gsub("MatureF(\\d|\\d{2})", "MatureF", pr$site)
pr$site <- gsub("ClimaxF(\\d|\\d{2})", "ClimaxF", pr$site)

pr$site <- as.factor(pr$site)
pr <- mutate(pr, site = factor(
  site, levels = c("DuneTop", "LeeSlo","InterLow", "WindSlo", "YoungF", 
                   "MatureF", "ClimaxF")))
ggthemr('solarized')
ggthemr_reset()
pr_plt <- 
  ggplot(data = pr, aes(x = site, y = projection)) +
  #geom_violin(scale = "count", aes(fill = site), 
  #           size = 0.3, width = 1.8) +
  geom_boxplot(width = 0.4, fill = NA, size = 0.3, 
               weight = 1) +
  labs(fill="site", shape="site", colour="site") +
  labs(x = "site", y = "Coverage") +
  ggtitle("Projective coverage grass") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

plot(pr_plt)
#stat_summary(fun = "median", geom = "point")
ggsave(filename = "prG_plt_v4.png", plot = pr_plt, dpi = 600, units = "in", 
       width = 7, height = 5)

#Lichens
dfL <- read.csv("Coverage - Lichn.csv", row.names = 1)

projection <- dfL[,1]
pr <- dfL

pr$site <- row.names(pr)
pr$site <- gsub("DuneTop(\\d|\\d{2})", "DuneTop", pr$site)
pr$site <- gsub("LeeSlo(\\d|\\d{2})", "LeeSlo", pr$site)
pr$site <- gsub("InterLow(\\d|\\d{2})", "InterLow", pr$site)
pr$site <- gsub("WindSlo(\\d|\\d{2})", "WindSlo", pr$site)
pr$site <- gsub("YoungF(\\d|\\d{2})", "YoungF", pr$site)
pr$site <- gsub("MatureF(\\d|\\d{2})", "MatureF", pr$site)
pr$site <- gsub("ClimaxF(\\d|\\d{2})", "ClimaxF", pr$site)

pr$site <- as.factor(pr$site)
pr <- mutate(pr, site = factor(
  site, levels = c("DuneTop", "LeeSlo","InterLow", "WindSlo", "YoungF", 
                   "MatureF", "ClimaxF")))
ggthemr('solarized')
ggthemr_reset()
pr_plt <- 
  ggplot(data = pr, aes(x = site, y = projection)) +
  #geom_violin(scale = "count", aes(fill = site), 
  #           size = 0.3, width = 1.8) +
  geom_boxplot(width = 0.4, fill = NA, size = 0.3, 
               weight = 1) +
  labs(fill="site", shape="site", colour="site") +
  labs(x = "site", y = "Coverage") +
  ggtitle("Projective coverage lichens") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

plot(pr_plt)
#stat_summary(fun = "median", geom = "point")
ggsave(filename = "prL_plt_v4.png", plot = pr_plt, dpi = 600, units = "in", 
       width = 7, height = 5)

#Mosses
dfM <- read.csv("Coverage - Moss.csv", row.names = 1)

projection <- dfM[,1]
pr <- dfM

pr$site <- row.names(pr)
pr$site <- gsub("DuneTop(\\d|\\d{2})", "DuneTop", pr$site)
pr$site <- gsub("LeeSlo(\\d|\\d{2})", "LeeSlo", pr$site)
pr$site <- gsub("InterLow(\\d|\\d{2})", "InterLow", pr$site)
pr$site <- gsub("WindSlo(\\d|\\d{2})", "WindSlo", pr$site)
pr$site <- gsub("YoungF(\\d|\\d{2})", "YoungF", pr$site)
pr$site <- gsub("MatureF(\\d|\\d{2})", "MatureF", pr$site)
pr$site <- gsub("ClimaxF(\\d|\\d{2})", "ClimaxF", pr$site)

pr$site <- as.factor(pr$site)
pr <- mutate(pr, site = factor(
  site, levels = c("DuneTop", "LeeSlo","InterLow", "WindSlo", "YoungF", 
                   "MatureF", "ClimaxF")))
ggthemr('solarized')
ggthemr_reset()
pr_plt <- 
  ggplot(data = pr, aes(x = site, y = projection)) +
  #geom_violin(scale = "count", aes(fill = site), 
  #           size = 0.3, width = 1.8) +
  geom_boxplot(width = 0.4, fill = NA, size = 0.3, 
               weight = 1) +
  labs(fill="site", shape="site", colour="site") +
  labs(x = "site", y = "Coverage") +
  ggtitle("Projective coverage mosses") +
  theme(legend.position = "none", plot.title = element_text(size = 12))

plot(pr_plt)
#stat_summary(fun = "median", geom = "point")
ggsave(filename = "prM_plt_v4.png", plot = pr_plt, dpi = 600, units = "in", 
       width = 7, height = 5)
