### Step 2 of TRIM-IT: K-MEANS CLUSTERING

##...............Function definition ....................
#Run the functions below before executing the main code:
# Data normalization
log_norm <-  function(x) log2(x+1)
#UMAP plot
plot.UMAP <- function(x, labels,
                      main="",
                      colors = custom,
                      pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
                      cex.main=1, cex.legend=0.85) {
  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  }
  
  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u <- unique(labels)
  legend.pos <- "topright"
  legend.text <- as.character(labels.u)
  if (add) {
    legend.pos <- "bottomright"
    legend.text <- paste(as.character(labels.u), legend.suffix)
  }
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

## .................Main Code.....................
# Library section
library(readr)
library("factoextra") #find the optimal number of cluster K
library(stats) #K-means clustering
library(cluster) #silhouette
library("fpc") #C-H score 
library(umap) 

# Load data:
load("~/Data/RNAseq-protein_informed-2021WHO.RData")
full.DS=GBM_RNA[order(rownames(GBM_RNA)),]
# Data normalization
full.DS=apply(full.DS, 2, log_norm)

# Find the optimal K:
set.seed(123)
fviz_nbclust(full.DS, kmeans, method='silhouette') # average silhouette 
#set.seed(123)
#fviz_nbclust(full.DS, kmeans, method='wss',print.summary=T) #+
 # geom_vline(xintercept = 3, linetype = 2)
set.seed(123)
fviz_nbclust(full.DS, kmeans, method='gap_stat') #gap statistic 


# K-means clustering applied to the GBM dataset
k=3 #optimal number of clusters
set.seed(123)
res.km=kmeans(full.DS, k, iter.max = 200, nstart = 20) #K-means

#Evaluate clustering performances
#silhouette:
d_t <- full.DS
dis <- dist(d_t)
si <- as.data.frame(silhouette(as.numeric(res.km$cluster), dis))
mean(si[,3])
#Calinski-Harabasz:
calinhara(full.DS,as.numeric(res.km$cluster)) 

# K-means clustering applied to the MOMIP-reduced GBM dataset
#load MOMIP results
load("~/Results/MOMIP-GBM-results.RData")
#reduce the dataset for clustering
MOMIP.DS=full.DS[,which(colnames(full.DS) %in% var_sel.2)]

set.seed(123)
k=3
res.km.MOMIP=kmeans(MOMIP.DS, k, iter.max = 200, nstart = 20) #K-means
#save(res.km.MOMIP,file="Kmeans-GBM-results.RData") #uncomment this line to save your results

#Evaluate clustering performances
#silhouette:
d_t <- MOMIP.DS
dis <- dist(d_t)
set.seed(15)
app <- silhouette(as.numeric(res.km.MOMIP$cluster), dis)
mean(app[,3])
#Calinski-Harabasz:
calinhara(MOMIP.DS,as.numeric(res.km.MOMIP$cluster)) 

#Cross-comparison between clusters 
table(res.km$cluster,res.km.MOMIP$cluster)

#UMAP:
set.seed(123)
full.DS.umap<- umap(full.DS)
MOMIP.umap<- umap(MOMIP.DS)

#plotting 
custom_colors <- c("C1" = "#FF6247" , "C2" = "#3CB371", "C3" = "#1E90FF")
plot.UMAP(full.DS.umap, as.numeric(res.km$cluster), colors = custom_colors)
title(main="full-DS")
plot.UMAP(MOMIP.umap, as.numeric(res.km.MOMIP$cluster), colors = custom_colors)
title(main="MOMIP-clusters")

