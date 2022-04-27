library(vegan)
library(tidyr)
library(dplyr)
library(R.utils)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(ggalt)
library(asbio)
library(ggthemr)
library(FSA)
library(rcompanion)
library(forcats)

setwd("/Users/admin/Yandex.Disk.localized/AnalisysSoromotin/") # setting the 
# working dirrectory

df <- read.csv("GeoBotany_fic.sp_newName_newLocacion.csv", row.names = 1)
# loading a dataframe with initial data

write.csv(df, 'GeoBotany_fic.sp_newName.csv')# writing dataframe to text files

# df$Fictitious.species <-1 #adding a dummy species

# doing a hierarchical cluster analysis of geobotanical data
dist_bray <- vegdist(df) # the function computes dissimilarity Bray-Curtis indices
hc <- hclust(dist_bray, method = "ward.D2") # hierarchical cluster analysis 
hc_geo_bot_plt <-   # on a set of dissimilarities and methods for analyzing it.
  fviz_dend(hc, k = 7,                 # Cut in six groups
            cex = 0.5,                 # label size
            color_labels_by_k = TRUE,  # color labels by groups
            ggtheme = theme_bw()     # Change theme
  )
plot(hc_geo_bot_plt) # displaying a graph based on cluster analysis
ggsave(filename = "hc_geo_bot_nL.png", plot = hc_geo_bot_plt, dpi = 600, 
       units = "in", width = 7, height = 5) # saving a graph to a file

##nMDS plotting based on geobotanical data
nMDS <- metaMDS(df, trymax = 100, distance = "bray", autotransform = FALSE)
# this function performs Nonmetric Multidimensional Scaling (NMDS), and tries 
# to find a stable solution using several random starts.

data.scores <- as.data.frame(scores(nMDS))
data.scores$site <- rownames(data.scores)
data.scores$grp <- rownames(data.scores)
data.scores$grp <- gsub("DT(\\d|\\d{2})", "DT", data.scores$grp)
data.scores$grp <- gsub("LS(\\d|\\d{2})", "LS", data.scores$grp)
data.scores$grp <- gsub("DL(\\d|\\d{2})", "DL", data.scores$grp)
data.scores$grp <- gsub("WS(\\d|\\d{2})", "WS", data.scores$grp)
data.scores$grp <- gsub("YF(\\d|\\d{2})", "YF", data.scores$grp)
data.scores$grp <- gsub("FF(\\d|\\d{2})", "FF", data.scores$grp)
data.scores$grp <- gsub("FG(\\d|\\d{2})", "FG", data.scores$grp)# creating 
# a grouping variable


species.scores <- as.data.frame(scores(nMDS, "species"))
species.scores$species <- rownames(species.scores)
grp.DT <- data.scores[data.scores$grp == "DT", ][chull(data.scores[
  data.scores$grp == "DT", c("NMDS1", "NMDS2")]), ]  # hull values for grp DT
grp.LS <- data.scores[data.scores$grp == "LS", ][chull(data.scores[
  data.scores$grp == "LS", c("NMDS1", "NMDS2")]), ]  # hull values for grp LS
grp.DL <- data.scores[data.scores$grp == "DL", ][chull(data.scores[
  data.scores$grp == "DL", c("NMDS1", "NMDS2")]), ]  # hull values for grp DL
grp.WS <- data.scores[data.scores$grp == "WS", ][chull(data.scores[
  data.scores$grp == "WS", c("NMDS1", "NMDS2")]), ]  # hull values for grp WS
grp.YF <- data.scores[data.scores$grp == "YF", ][chull(data.scores[
  data.scores$grp == "YF", c("NMDS1", "NMDS2")]), ]  # hull values for grp YF
grp.FF <- data.scores[data.scores$grp == "FF", ][chull(data.scores[
  data.scores$grp == "FF", c("NMDS1", "NMDS2")]), ]  # hull values for grp FF
grp.FG <- data.scores[data.scores$grp == "FG", ][chull(data.scores[
  data.scores$grp == "FG", c("NMDS1", "NMDS2")]), ]  # hull values for grp FG
hull.data <- rbind(grp.DT, grp.LS, grp.DL, grp.WS, grp.YF, grp.FF, grp.FG)#combine 
# groups

geo_bot_plt <-
  ggplot() + # creating an nMDS graph
  geom_encircle(data = data.scores, aes(x=NMDS1,y=NMDS2, fill = grp), 
                alpha=0.4) +   # draw circles
  #geom_label_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), 
                   #colour = "red", box.padding = unit(0.04, "lines"), force = 4, 
                   #max.overlaps = Inf, segment.color = NA,
                   #size=2.5, alpha=0.8) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp, colour=grp),
             size=1.5) + # add the point markers
  scale_shape_manual(values=(1:7)) + # change shape the point markers
  labs(fill="Sites")+
  labs(shape="Sites", colour="Sites")+
  guides(fill = guide_legend(override.aes = list(linetype = 0))) +
  geom_text_repel(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),
  box.padding = unit(0.04, "lines"), force_pull = 3, force = 4, 
  min.segment.lengt=0.4, segment.size = 0.3, segment.alpha = 0.7,
  nudge_x=0.06, nudge_y=0.04, size=2.5, 
  max.overlaps = Inf) +  # add the site labels
  coord_equal() +
  theme_bw()+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=10), # remove x-axis labels
        axis.title.y = element_text(size=10), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
plot(geo_bot_plt)

ggsave(filename = "nMDS_geo_bot(circle)_3var.png", plot = geo_bot_plt, 
       bg='#ffffff', dpi = 600, 
       units = "in", width = 7, height = 5) # saving a graph to a file

#PERMANOVA

#create a grouping variable in the main dataframe and assign it as a factor.
df$Group <- row.names(df)
df$Group <- gsub("\\d", "", df$Group)
df$Group <- as.factor(df$Group)

# Since the groups have the same volume, the equality of intra-group variances 
# can not be checked 

# permanova procedure
permanova <- adonis(df[, 1:ncol(df)-1] ~ Group, data = df, 
  permutations = 999, method = "bray")
permanova
