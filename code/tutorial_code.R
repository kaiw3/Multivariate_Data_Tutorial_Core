library(vegan)
library(dplyr)
library(gplots)
library(tidyverse)

barents <- read.csv("data/barents_data.csv")
barents_spp <- barents %>% select(Re_hi:Tr_spp)
barents_env_raw <- barents %>% select(ID.No:Temperature)

barents_env <- barents_env_raw %>% mutate(Depth_cat=cut(Depth, breaks=c(-Inf, 250, 350, Inf), labels=c("shallow","mid","deep"))) %>%
  mutate(Temp_cat=cut(Temperature, breaks=c(-Inf, 1, 2, 3, Inf), labels=c("vc","c","m","w")))


barents.mds <- metaMDS(barents_spp, distance = "bray", autotransform = FALSE)
barents.envfit <- envfit(barents.mds, barents_env, permutations = 999) # this fits environmental vectors
barents.sppfit <- envfit(barents.mds, barents_spp, permutations = 999) # this fits species vectors
barents.mds  # Stress value is less than 0.2, which is good. Shows how easy it was to condense multidimensional data into two dimensional space

site.scrs <- as.data.frame(scores(barents.mds, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, Depth = barents_env$Depth_cat) #add grouping variable "Management" to dataframe
site.scrs <- cbind(site.scrs, Temperature = barents_env$Temp_cat) #add grouping variable of cluster grouping to dataframe
# site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot
head(site.scrs)

# USED TO SHOW SPECIES DATA ON PLOT IF YOU WANT
spp.scrs <- as.data.frame(scores(barents.sppfit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = barents.sppfit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05
head(spp.scrs)

# Environmental extrinsic variables
env.scores.barents <- as.data.frame(scores(barents.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores.barents <- cbind(env.scores.barents, env.variables = rownames(env.scores.barents)) #and then gives them their names

env.scores.barents <- cbind(env.scores.barents, pval = barents.envfit$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores.barents, pval<=0.05) #subset data to show variables significant at 0.05

head(env.scores.barents)

# Basic nmds
nmds.plot.barents <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(x = NMDS1, y = NMDS2, colour = factor(site.scrs$Depth), shape = factor(site.scrs$Temperature)), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Depth", shape = "Temperature")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

nmds.plot.barents + labs(title = "Basic ordination plot") #displays plot


nmds.plot.barents <- nmds.plot.barents+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with species vectors")

nmds.plot.barents+
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)+ #add labels for env variables
  labs(title="Ordination with environmental vectors")

# We can see that theres a distinction in samples based on depth and temperature. 
# Species info shows some species are more commonly found in different conditions.
# None of this is statistically verified yet, so a ANOSIM, SIMPER or indicator species analysis will need to be run to check.








barents_pca$x
barents_pca <- prcomp(barents_env_raw)
summary(barents_pca)
screeplot(barents_pca, type = "lines")
(barents_pca$sdev)^2

y <- matrix(rnorm(89)) 
c <- cor(barents_spp, method="spearman")
d <- as.dist(1-c)
pear.dist <- as.dist((1-c)/2)
levelplot(c)
hp <- hclust(pear.dist, method="ward.D")
plot(hp, hang=-1, cex=0.2)

data.med.center <- sweep(barents_spp, MARGIN=1, STATS=apply(barents_spp,1,median), FUN = "-",) 
data.med.center[ data.med.center > 5 ] <- 5
data.med.center[ data.med.center < -5 ] <- -5

pear.cols <- cor(t(data.med.center), method="pearson")
pear.dist <- as.dist((1-pear.cols)/2)
hc.cols <- hclust(pear.dist, method="ward.D")
pear.rows <- cor(data.med.center, method="pearson")
pear.dist <- as.dist((1-pear.rows)/2)
hc.rows <- hclust(pear.dist, method="ward.D")

heatmap.2(t(data.med.center), 
          scale="none", 
          cluster.by.col=TRUE, 
          cluster.by.row=TRUE, 
          Rowv=TRUE,
          Colv=TRUE,
          hclust.col = hc.cols, 
          hclust.row = hc.rows, 
          labCol=NULL,
          cexRow = 0.5,
          )


# Failed NMDS
community_matrix <- matrix(
  sample(1:100,300,replace=T),nrow=10,
  dimnames=list(paste("community",1:10,sep=""),paste("sp",1:30,sep="")))

example_NMDS <- metaMDS(barents_spp, k = 2)
stressplot(example_NMDS)        # Produce shepherds diagram to test goodness of fit

plot(example_NMDS, "sites")   # Produces distance 
orditorp(example_NMDS, "sites")
colvec <- c("gray0", "gray0", "gray49", "gray49")   # Identifies colors for group assignments
pchvec <- c(21, 21, 22, 22)

ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)







