# First attempts at selecting a test to use in the tutorial

# Load packages
library(vegan)  # For performing the nmds
library(gplots)  # For heatmaps
library(tidyverse)  # For data tidying functions and plotting

# Load data
barents <- read.csv("data/barents_data.csv")

barents_pca <- prcomp(barents_env_raw)
barents_pca$x
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



