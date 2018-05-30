#20180530_pull_indiv_clusters.R
#Olga Minkina
#Pull only sme clusters from k=10, k=20 clustered data.

#load Barcode-Gene-Cluster Color_DataFrame for clusters 2-10 and graph-based clusters
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180522 Analysis_scripts and outputs/20180522_script_2_results/")
df_10 <- read.csv("Barcode-Gene-Cluster-Color_DataFrame_20180522.csv")

#load Barcode-Gene-Cluster-Color_DataFrame for clusters 11-20
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180530_k1120Clusters/Resulting_cluster_plots/")
df_20 <- read.csv("Barcode-Gene-Cluster-Color_DataFrame_20180530.csv")

#Get only relevant columns 
df_10_clust10only <- df_10[, c(2, 3, 4, 23, 24, 34)]
df_20_clust20only <- df_20[, c(2, 3, 4, 23, 24, 34)]

#plot k=10, k=20 clusters (just to make sure df is correct)
plot(df_10_clust10only$TSNE.1, df_10_clust10only$TSNE.2, col = paste(df_10_clust10only[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("k=10", adj = 0.05, line = -2))
plot(df_20_clust20only$TSNE.1, df_20_clust20only$TSNE.2, col = paste(df_20_clust20only[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("k=20", adj = 0.05, line = -2))

#get only rows corresponding to clusters of interest
#plot to make sure I got the correct cluster
#k=10
condition <- df_10_clust10only$kmeans_10_Cluster
k10_yto<-df_10_clust10only[condition == 2 | condition ==3 | condition ==5,]
plot(k10_yto$TSNE.1, k10_yto$TSNE.2, col = paste(k10_yto[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("k=10, pulled clusters", adj = 0.05, line = -2), ylim=c(-50,50), xlim = c(-50, 50))
dev.copy(png, 'cluster_kmeans_10_clusterspulled.png', width = 500, height = 500)
dev.off()
#k=20
condition <- df_20_clust20only$kmeans_20_Cluster
k20_yto<-df_20_clust20only[condition == 1 | condition ==8 | condition ==6 | condition ==12 |condition ==5,]
plot(k20_yto$TSNE.1, k20_yto$TSNE.2, col = paste(k20_yto[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("k=20, pulled clusters", adj = 0.05, line = -2), ylim=c(-50,50), xlim = c(-50, 50))
dev.copy(png, 'cluster_kmeans_20_clusterspulled.png', width = 500, height = 500)
dev.off()

