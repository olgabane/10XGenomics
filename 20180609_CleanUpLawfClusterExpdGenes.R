#20180609_CleanUpLawfClusterExpdGenes.R
#Delete genes with expression only in a few cells. 

#Load GBM file created 20180531 (see 20180531_GBMpulled.R)
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180531_GBMpulled/")
k20_Lawf_cluster_GBM_wGeneNames<-read.table("k20_Lawf_cluster_GBM_wGeneNames.txt")

###Delete genes that are not expressed in any of the 852 cells
#sum across rows
sum<-c()
for (i in 1:dim(k20_Lawf_cluster_GBM_wGeneNames)[1])
  sum[i] <- sum(k20_Lawf_cluster_GBM_wGeneNames[i, 4:855])
#add sums in column to dataframe
k20_Lawf_cluster_GBM_wGeneNames$RowSum <- sum
length(which(k20_Lawf_cluster_GBM_wGeneNames$RowSum == 0)) #3019 genes eliminated (not expressed in any cells in this cluster)
#delete all rows for which RowSum == 0. I did some manual check to make sure this worked properly.
zerorows<-which(k20_Lawf_cluster_GBM_wGeneNames$RowSum == 0)
k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly <- k20_Lawf_cluster_GBM_wGeneNames[-zerorows,]
dim(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly)
#Sort by expression (RowSums), in descending order
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180609_CleanUpLawfCluterExpdGenes/")
k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted <-
  k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[order(-k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly$RowSum),]
write.table(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted, 
            "k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted.txt", sep = "\t")

#plot: RowSums across 852 cells
quartz("gene expr sums", 6, 4.5)
plot(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$RowSum, 
     xlab = "8935 Genes", ylab = "Log10(Normalized UMI) for 852 cells", pch = 16, cex.lab = 0.7, 
     cex.axis = 0.7, cex = 0.5)



#####for each gene, now many cells have expression > 0?
