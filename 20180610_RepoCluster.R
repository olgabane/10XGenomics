#20180610_RepoCluster.R
#Processed repo+ cluster cells 
#Compare repo expressed cluster to Lawf

#Get repo+ data frame
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180531_GBMpulled/")
k20_repo_cluster_GBM_wGeneNames <- read.table("k20_repo_cluster_GBM_wGeneNames.txt")

###Delete genes that are not expressed in any of the 1122 cells
#sum across rows
sum<-c()
for (i in 1:dim(k20_repo_cluster_GBM_wGeneNames)[1])
  sum[i] <- sum(k20_repo_cluster_GBM_wGeneNames[i, 4:1122])
#add sums in column to dataframe
k20_repo_cluster_GBM_wGeneNames$RowSum <- sum
length(which(k20_repo_cluster_GBM_wGeneNames$RowSum == 0)) #2130 genes eliminated (not expressed in any cells in this cluster)
#delete all rows for which RowSum == 0. I did some manual check to make sure this worked properly.
zerorows<-which(k20_repo_cluster_GBM_wGeneNames$RowSum == 0)
k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly <- k20_repo_cluster_GBM_wGeneNames[-zerorows,]
#Sort by expression (RowSums), in descending order
k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted <-
  k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly[order(-k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly$RowSum),]

###plot: RowSums across 1122 cells
quartz("gene expr sums", 6, 4.5)
plot(k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$RowSum, 
     xlab = "1122 Genes", ylab = "Log10(Normalized UMI) for 1122 cells", pch = 16, cex.lab = 0.7, 
     cex.axis = 0.7, cex = 0.5)

###for each gene, now many cells have expression > 0?
zerocount <-c()
expressedcount <-c()
for (i in 1:dim(k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted)[1])
{
  a<-as.vector(k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted[i, 4:1122])
  zerocount[i] <- sum(a == 0) 
  expressedcount[i] <- sum(a > 0)
}

k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$Cells0 <- zerocount
k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0 <- expressedcount

quartz("gene expr sums", 6, 4.5)
plot(k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0, 
     xlab = "8935 Genes", ylab = "Number of cells that express gene (>0)", pch = 16, cex.lab = 0.7, 
     cex.axis = 0.7, cex = 0.5)
abline(h=50, col = "red")

cutoff <- tail(which(k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0 >= 50), 1)

#remove genes with <50 cells expressing gene
k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted_top6272 <- 
  k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted[-c(cutoff+1:9824), ]

top6272genes <- as.vector(k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted_top6272$FBN_Flybase)
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180610_repolawfclustercompared/")
write(top6272genes, "top6272genes_repo_cluster.txt")

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180609_CleanUpLawfCluterExpdGenes/")
top5304genes <- read.table("top5304genes.txt")

#compare 6272 genes exp'd in repo cluster with 5304 genes expressed in L cluster.
shared_genes <- intersect(top6272genes, top5304genes$V1)



