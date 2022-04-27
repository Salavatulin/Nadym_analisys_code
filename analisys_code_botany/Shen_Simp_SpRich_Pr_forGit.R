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
simp$site <- gsub("DT(\\d|\\d{2})", "DT", simp$site)
simp$site <- gsub("LS(\\d|\\d{2})", "LS", simp$site)
simp$site <- gsub("DL(\\d|\\d{2})", "DL", simp$site)
simp$site <- gsub("WS(\\d|\\d{2})", "WS", simp$site)
simp$site <- gsub("YF(\\d|\\d{2})", "YF", simp$site)
simp$site <- gsub("FF(\\d|\\d{2})", "FF", simp$site)
simp$site <- gsub("FG(\\d|\\d{2})", "FG", simp$site)

simp$site <- as.factor(simp$site)
simp <- mutate(simp, site = factor(
  site, levels = c("DT", "LS","DL", "WS", "YF", "FF", "FG")))
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
ggsave(filename = "simp_plt_2.png", plot = simp_plt, dpi = 600, units = "in", 
       width = 7, height = 5)

# Shannon
shannon <- diversity(df, index = "shannon")
shan <- as.data.frame(shannon)
shan$site <- row.names(shan)
shan$site <- gsub("DT(\\d|\\d{2})", "DT", shan$site)
shan$site <- gsub("LS(\\d|\\d{2})", "LS", shan$site)
shan$site <- gsub("DL(\\d|\\d{2})", "DL", shan$site)
shan$site <- gsub("WS(\\d|\\d{2})", "WS", shan$site)
shan$site <- gsub("YF(\\d|\\d{2})", "YF", shan$site)
shan$site <- gsub("FF(\\d|\\d{2})", "FF", shan$site)
shan$site <- gsub("FG(\\d|\\d{2})", "FG", shan$site)

shan$site <- as.factor(shan$site)
shan <- mutate(shan, site = factor(
  site, levels = c("DT", "LS","DL", "WS", "YF", "FF", "FG")))
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
ggsave(filename = "shan_plt_2.png", plot = shan_plt, dpi = 600, units = "in", 
       width = 7, height = 5)


#Species richness

richness <- specnumber(df)
rich <- as.data.frame(richness)

rich$site <- row.names(rich)
rich$site <- gsub("DT(\\d|\\d{2})", "DT", rich$site)
rich$site <- gsub("LS(\\d|\\d{2})", "LS", rich$site)
rich$site <- gsub("DL(\\d|\\d{2})", "DL", rich$site)
rich$site <- gsub("WS(\\d|\\d{2})", "WS", rich$site)
rich$site <- gsub("YF(\\d|\\d{2})", "YF", rich$site)
rich$site <- gsub("FF(\\d|\\d{2})", "FF", rich$site)
rich$site <- gsub("FG(\\d|\\d{2})", "FG", rich$site)

rich$site <- as.factor(rich$site)
rich <- mutate(rich, site = factor(
  site, levels = c("DT", "LS","DL", "WS", "YF", "FF", "FG")))
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
ggsave(filename = "rich_plt.png", plot = rich_plt, dpi = 600, units = "in", 
       width = 7, height = 5)

##Projective coverage

#Grass
dfG <- read.csv("Coverage - Grass.csv", row.names = 1)

projection <- dfG[,1]
pr <- dfG

pr$site <- row.names(pr)
pr$site <- gsub("DT(\\d|\\d{2})", "DT", pr$site)
pr$site <- gsub("LS(\\d|\\d{2})", "LS", pr$site)
pr$site <- gsub("DL(\\d|\\d{2})", "DL", pr$site)
pr$site <- gsub("WS(\\d|\\d{2})", "WS", pr$site)
pr$site <- gsub("YF(\\d|\\d{2})", "YF", pr$site)
pr$site <- gsub("FF(\\d|\\d{2})", "FF", pr$site)
pr$site <- gsub("FG(\\d|\\d{2})", "FG", pr$site)

pr$site <- as.factor(pr$site)
pr <- mutate(pr, site = factor(
  site, levels = c("DT", "LS","DL", "WS", "YF", "FF", "FG")))
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
ggsave(filename = "prG_plt.png", plot = pr_plt, dpi = 600, units = "in", 
       width = 7, height = 5)

#Lichens
dfL <- read.csv("Coverage - Lichn.csv", row.names = 1)

projection <- dfL[,1]
pr <- dfL

pr$site <- row.names(pr)
pr$site <- gsub("DT(\\d|\\d{2})", "DT", pr$site)
pr$site <- gsub("LS(\\d|\\d{2})", "LS", pr$site)
pr$site <- gsub("DL(\\d|\\d{2})", "DL", pr$site)
pr$site <- gsub("WS(\\d|\\d{2})", "WS", pr$site)
pr$site <- gsub("YF(\\d|\\d{2})", "YF", pr$site)
pr$site <- gsub("FF(\\d|\\d{2})", "FF", pr$site)
pr$site <- gsub("FG(\\d|\\d{2})", "FG", pr$site)

pr$site <- as.factor(pr$site)
pr <- mutate(pr, site = factor(
  site, levels = c("DT", "LS","DL", "WS", "YF", "FF", "FG")))
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
ggsave(filename = "prL_plt.png", plot = pr_plt, dpi = 600, units = "in", 
       width = 7, height = 5)

#Mosses
dfM <- read.csv("Coverage - Moss.csv", row.names = 1)

projection <- dfM[,1]
pr <- dfM

pr$site <- row.names(pr)
pr$site <- gsub("DT(\\d|\\d{2})", "DT", pr$site)
pr$site <- gsub("LS(\\d|\\d{2})", "LS", pr$site)
pr$site <- gsub("DL(\\d|\\d{2})", "DL", pr$site)
pr$site <- gsub("WS(\\d|\\d{2})", "WS", pr$site)
pr$site <- gsub("YF(\\d|\\d{2})", "YF", pr$site)
pr$site <- gsub("FF(\\d|\\d{2})", "FF", pr$site)
pr$site <- gsub("FG(\\d|\\d{2})", "FG", pr$site)

pr$site <- as.factor(pr$site)
pr <- mutate(pr, site = factor(
  site, levels = c("DT", "LS","DL", "WS", "YF", "FF", "FG")))
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
ggsave(filename = "prM_plt.png", plot = pr_plt, dpi = 600, units = "in", 
       width = 7, height = 5)
