#20180530_clusers.R
#Olga Minkina 
#Plot clusters (manually, withouth cellranger). Similar to 20180522 script.

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180528_Prince/Files from cellranger_reanalyze/outs/analysis/tsne/2_components/")
tsne_file <- read.csv("projection.csv")

#get barcodes and tsne components
for(i in 1:3)
  assign(paste0("x", i), tsne_file[i])
#get barcodes and clusters for k-means clustering 11-20
for(i in 11:20)
{
  setwd(paste0("../../clustering/kmeans_", i,"_clusters" ))
  cluster_file<-read.csv("clusters.csv")
  assign(paste0("x", i), cluster_file[,1:2])
}

#create dataframe with barcodes, TSNE components and cluster information
#nested loop could make this more streamlined
df = data.frame(Barcode = x1, TSNE.1 = x2, TSNE.2 = x3, 
                kmeans_11_Barcode = x11[, 1], kmeans_11_Cluster = x11[,2],
                kmeans_12_Barcode = x12[, 1], kmeans_12_Cluster = x12[,2],
                kmeans_13_Barcode = x13[, 1], kmeans_13_Cluster = x13[,2],
                kmeans_14_Barcode = x14[, 1], kmeans_14_Cluster = x14[,2],
                kmeans_15_Barcode = x15[, 1], kmeans_15_Cluster = x15[,2],
                kmeans_16_Barcode = x16[, 1], kmeans_16_Cluster = x16[,2],
                kmeans_17_Barcode = x17[, 1], kmeans_17_Cluster = x17[,2],
                kmeans_18_Barcode = x18[, 1], kmeans_18_Cluster = x18[,2],
                kmeans_19_Barcode = x19[, 1], kmeans_19_Cluster = x19[,2],
                kmeans_20_Barcode = x20[, 1], kmeans_20_Cluster = x20[,2])

#make vector of 20 colors
color_vect <- c("red", "darkorange", "yellow", "lightgreen", "tan3", "cyan", "dodgerblue", 
                "darkred", "darkorchid2", "gray53", "black", "pink2", "darkmagenta", "darkcyan", 
                "navy", "darkgreen", "darkseagreen1", "darkslateblue", "khaki", "yellow4") 

#make column names for color assignments
color_column_list <-list()
for(j in 11:20)
  color_column_list[j-10] = paste("kmeans_", j, "_Color", sep = "")

#add columns to dataframe to assign colors to clusters using k-means 11-20 clustering
for(j in 11:20)
  df[,color_column_list[[j-10]]] <- NA

#cluster column list
cluster_column_list <-list()
for(j in 11:20)
  cluster_column_list[j-10] = paste("kmeans_", j, "_Cluster", sep = "")

#assign colors to 2-10 clusters generated  using k-means 2-10 clustering 
for(i in 1:20)
  for(j in 2:11)
    df[,color_column_list[[j-1]]][df[,cluster_column_list[[j-1]]] == i] = color_vect[i]

#change working directory and save resulting dataframe
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180530_k1120Clusters/Resulting_cluster_plots/")
write.csv(df, "Barcode-Gene-Cluster-Color_DataFrame_20180530.csv")

#plot and save clusters generated  using k-means 11-20 clustering 
for(j in 2:11)
{
  png(paste("clusters_kmeans",j+9, ".png"))
  plot(df$TSNE.1, df$TSNE.2, col = df[,color_column_list[[j-1]]], xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title(paste("k =",j+9), adj = 0.05, line = -2))
  dev.off()
}

#compare these clusters to clusters plotted by cellRanger (should be exactly the same)
library(cellrangerRkit)
cellranger_pipestance_path <- ("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180528_Prince/Files from cellranger_reanalyze/")
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)
n_clu <- 11:20
km_res <- analysis_results$clustering # load pre-computed kmeans results
clu_res <- sapply(n_clu, function(x) km_res[[paste("kmeans",x,"clusters",sep="_")]]$Cluster)
colnames(clu_res) <- sapply(n_clu, function(x) paste("kmeans",x,sep="."))
tsne_proj <- analysis_results$tsne
png('kmeans_11to20_clusters_cellRanger.png', width = 800, height = 600)
visualize_clusters(clu_res,tsne_proj[c("TSNE.1","TSNE.2")])
dev.off()

#Save relevant plots in Demo folder
setwd("/Users/Olga/Desktop/Demo/20180530")
#plot and save clusters generated  using k-means 11-20 clustering 
for(j in 2:11)
{
  png(paste("clusters_kmeans",j+9, ".png"))
  plot(df$TSNE.1, df$TSNE.2, col = df[,color_column_list[[j-1]]], xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title(paste("k =",j+9), adj = 0.05, line = -2))
  dev.off()
}
