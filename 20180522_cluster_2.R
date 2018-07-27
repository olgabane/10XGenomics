#20180522_script_2.R
#Olga Minkina 
#May 22, 2018
#Does same thing as 20180522_cluster_1.R, but gets data in less obfuscated way (without using cell ranger). Also more stream-lined.

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/cellranger/count-Neuron/outs/analysis/")
setwd("tsne/2_components/")
tsne_file <- read.csv("projection.csv")

#get barcodes and tsne components
for(i in 1:3)
  assign(paste0("x", i), tsne_file[i])
#get barcodes and clusters for graph based k-means clustering
setwd("../../clustering/graphclust/")
cluster_file<-read.csv("clusters.csv")
x4 <- cluster_file[,1:2]
#get barcodes and clusters for k-means clustering 2-10
for(i in 2:10)
{
  setwd(paste0("../../clustering/kmeans_", i,"_clusters" ))
  cluster_file<-read.csv("clusters.csv")
  assign(paste0("x", i+3), cluster_file[,1:2])
}

#create dataframe with barcodes, TSNE components and cluster information
#nested loop could make this more streamlined
df = data.frame(Barcode = x1, TSNE.1 = x2, TSNE.2 = x3, 
                graphclust_Barcode = x4[,1], graphclust_Cluster = x4[,2],
                kmeans_2_Barcode = x5[, 1], kmeans_2_Cluster = x5[,2],
                kmeans_3_Barcode = x6[, 1], kmeans_3_Cluster = x6[,2],
                kmeans_4_Barcode = x7[, 1], kmeans_4_Cluster = x7[,2],
                kmeans_5_Barcode = x8[, 1], kmeans_5_Cluster = x8[,2],
                kmeans_6_Barcode = x9[, 1], kmeans_6_Cluster = x9[,2],
                kmeans_7_Barcode = x10[, 1], kmeans_7_Cluster = x10[,2],
                kmeans_8_Barcode = x11[, 1], kmeans_8_Cluster = x11[,2],
                kmeans_9_Barcode = x12[, 1], kmeans_9_Cluster = x12[,2],
                kmeans_10_Barcode = x13[, 1], kmeans_10_Cluster = x13[,2])

###Assign colors to clusters 1-16
#make vector of 16 colors
color_vect <- c("red", "darkorange", "yellow", "lightgreen", "tan3", "cyan", "dodgerblue", 
                "darkred", "darkorchid2", "gray53", "black", "pink2", "darkmagenta", "darkcyan", 
                "navy", "darkgreen")

#assign colors to up to 16 clusters made from graph-based k-means clustering
for(i in 1:16)
  df$graphclust_Color[df$graphclust_Cluster == i] = color_vect[i]

#make column names for color assignments for k= 2-10
color_column_list <-list()
for(j in 2:10)
  color_column_list[j-1] = paste("kmeans_", j, "_Color", sep = "")
#add columns to dataframe to assign colors to clusters using k-means 2-10 clustering
for(j in 2:10)
  df[,color_column_list[[j-1]]] <- NA

#cluster column name list
cluster_column_list <-list()
for(j in 2:10)
  cluster_column_list[j-1] = paste("kmeans_", j, "_Cluster", sep = "")

#assign colors to 2-10 clusters generated  using k-means 2-10 clustering 
for(i in 1:16)
  for(j in 2:10)
    df[,color_column_list[[j-1]]][df[,cluster_column_list[[j-1]]] == i] = color_vect[i]

#change working directory and save resulting dataframe
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180522 Analysis_scripts and outputs/20180522_script_2_results/")
write.csv(df, "Barcode-Gene-Cluster-Color_DataFrame_20180522.csv")

#plot and save graph-based clusters
png("clusters_kmeans_graph.png")
plot(df$TSNE.1, df$TSNE.2, col = df$graphclust_Color, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("Graph based", adj = 0.05, line =-2))
dev.off()

#plot and save clusters generated  using k-means 2-10 clustering 
for(j in 2:10)
{
  png(paste("clusters_kmeans",j, ".png"))
  plot(df$TSNE.1, df$TSNE.2, col = df[,color_column_list[[j-1]]], xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title(paste("k =",j), adj = 0.05, line = -2))
  dev.off()
}

#Save relevant plots in Demo folder
setwd("/Users/Olga/Desktop/Demo/20180522")
#plot and save clusters generated using graph-based clustering and k-means 2-10 clustering 
png("clusters_kmeans_graph.png")
plot(df$TSNE.1, df$TSNE.2, col = df$graphclust_Color, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("Graph based", adj = 0.05, line =-2))
dev.off()

for(j in 2:10)
{
  png(paste("clusters_kmeans",j, ".png"))
  plot(df$TSNE.1, df$TSNE.2, col = df[,color_column_list[[j-1]]], xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title(paste("k =",j), adj = 0.05, line = -2))
  dev.off()
}


