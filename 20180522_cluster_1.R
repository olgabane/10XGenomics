#20180522_script_1.R
#Olga Minkina 
#May 22, 2018
#Initial script to plot clusters on top of tsne

#load cell ranger and get analysis_results, a file that aggregates tsne and cluster information
library(cellrangerRkit)
cellranger_pipestance_path <- "/Users/Olga/Downloads/cellranger/count-Neuron"
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)

#build data frame with barcode and cluster information
#this might be worth putting in a few for loops for nicer code
x1 <-  analysis_results$tsne$Barcode #NOTE: same as "analysis_results$tsne[,"Barcode"]" and "analysis_results$tsne[,1]"
x2 <-  analysis_results$tsne$TSNE.1
x3 <-  analysis_results$tsne$TSNE.2
x4 <- analysis_results$clustering$graphclust[,1:2]
x5<-analysis_results$clustering$kmeans_10_clusters[,1:2]
x6<-analysis_results$clustering$kmeans_9_clusters[,1:2]
x7<-analysis_results$clustering$kmeans_8_clusters[,1:2]
x8<-analysis_results$clustering$kmeans_7_clusters[,1:2]
x9<-analysis_results$clustering$kmeans_6_clusters[,1:2]
x10<-analysis_results$clustering$kmeans_5_clusters[,1:2]
x11<-analysis_results$clustering$kmeans_4_clusters[,1:2]
x12<-analysis_results$clustering$kmeans_3_clusters[,1:2]
x13<-analysis_results$clustering$kmeans_2_clusters[,1:2]
df = data.frame(Barcode = x1, TSNE.1 = x2, TSNE.2 = x3, 
                graphclust_Barcode = x4[,1], graphclust_Cluster = x4[,2],
                kmeans_10_Barcode = x5[, 1], kmeans_10_Cluster = x5[,2],
                kmeans_9_Barcode = x6[, 1], kmeans_9_Cluster = x6[,2],
                kmeans_8_Barcode = x7[, 1], kmeans_8_Cluster = x7[,2],
                kmeans_7_Barcode = x8[, 1], kmeans_7_Cluster = x8[,2],
                kmeans_6_Barcode = x9[, 1], kmeans_6_Cluster = x9[,2],
                kmeans_5_Barcode = x10[, 1], kmeans_5_Cluster = x10[,2],
                kmeans_4_Barcode = x11[, 1], kmeans_4_Cluster = x11[,2],
                kmeans_3_Barcode = x12[, 1], kmeans_3_Cluster = x12[,2],
                kmeans_2_Barcode = x13[, 1], kmeans_2_Cluster = x13[,2])

#make column names for color assignments
color_column_list <-list()
for(j in 2:10)
  color_column_list[j-1] = paste("kmeans_", j, "_Color", sep = "")

#make vector of 16 colors
color_vect <- c("red", "darkorange", "yellow", "lightgreen", "tan3", "cyan", "dodgerblue", 
                  "darkred", "darkorchid2", "gray53", "black", "pink2", "darkmagenta", "darkcyan", 
                  "navy", "darkgreen")

#assign colors to 16 clusters made from graph-based k-means clustering
for(i in 1:16)
  df$graphclust_Color[df$graphclust_Cluster == i] = color_vect[i]

#add columns to dataframe to assign colors to clusters using k-means 2-10 clustering
for(j in 2:10)
  df[,color_column_list[[j-1]]] <- NA

#cluster column list
cluster_column_list <-list()
for(j in 2:10)
  cluster_column_list[j-1] = paste("kmeans_", j, "_Cluster", sep = "")

#assign colors to 2-10 clusters generated  using k-means 2-10 clustering 
for(i in 1:16)
{
  for(j in 2:10)
    df[,color_column_list[[j-1]]][df[,cluster_column_list[[j-1]]] == i] = color_vect[i]
}

#change working directory and save resulting dataframe
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180522 Analysis_scripts and outputs/20180522_script_1_results/")
write.csv(df, "Barcode-Gene-Cluster-Color_DataFrame_20180522.csv")

#plot and save graph-based clusters
png("clusters_kmeans_graph.png")
plot(df$TSNE.1, df$TSNE.2, col = df$graphclust_Color, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5)
dev.off()

#plot and save clusters generated  using k-means 2-10 clustering 
for(j in 2:10)
{
  filename <- paste("clusters_kmeans",j, ".png")
  png(filename)
  plot(df$TSNE.1, df$TSNE.2, col = df[,color_column_list[[j-1]]], xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5)
  dev.off()
}
