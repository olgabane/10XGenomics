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

###plot: RowSums across 852 cells
quartz("gene expr sums", 6, 4.5)
plot(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$RowSum, 
     xlab = "8935 Genes", ylab = "Log10(Normalized UMI) for 852 cells", pch = 16, cex.lab = 0.7, 
     cex.axis = 0.7, cex = 0.5)

which(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted[1, 4:855] == 0)

###for each gene, now many cells have expression > 0?
zerocount <-c()
expressedcount <-c()
for (i in 1:dim(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted)[1])
{
  a<-as.vector(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted[i, 4:855])
  zerocount[i] <- sum(a == 0) 
  expressedcount[i] <- sum(a > 0)
}

k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$Cells0 <- zerocount
k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0 <- expressedcount

quartz("gene expr sums", 6, 4.5)
plot(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0, 
     xlab = "8935 Genes", ylab = "Number of cells that express gene (>0)", pch = 16, cex.lab = 0.7, 
     cex.axis = 0.7, cex = 0.5)
abline(h=50, col = "red")

#plots both graphs on top of each other
png("SumGeneExp.png", width = 600, height = 500)
layout(matrix(c(1,2), nrow = 2, byrow = TRUE))
par(mar = c(4, 2, 1, 1), oma = c (0, 4, 1, 1), xpd = NA)
plot(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$RowSum,
     xlab = "8935 Genes" , ylab = "Summed Expression", pch = 16, cex = 0.5, cex.lab = 1)
plot(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0, 
     xlab = "8935 Genes", ylab = "Cells that express gene (>0)", pch = 16, cex = 0.5, cex.lab = 1)
#par(xpd = TRUE)
lines(x = c(0, 8935), y = c(50, 50), lty = 1, col = "red", lwd = 2)
text(250, 100, "y=50", col = "red")
#abline(h = 50, col = "red")
dev.off()


#re-plot highlighting genes I will keep
cutoff <- tail(which(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0 >= 50), 1)
png("SumGeneExp_colored.png", width = 600, height = 500)
layout(matrix(c(1,2), nrow = 2, byrow = TRUE))
par(mar = c(4, 2, 1, 1), oma = c (0, 4, 1, 1), xpd = NA)
plot(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$RowSum,
     xlab = "8935 Genes" , ylab = "Summed Expression", pch = 16, cex = 0.5, cex.lab = 1, 
     col = ifelse(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0 >= 50, "blue", "black"))
plot(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0, 
     xlab = "8935 Genes", ylab = "Cells that express gene (>0)", pch = 16, cex = 0.5, cex.lab = 1, 
     col = ifelse(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted$CellsGreaterThan0 >= 50, "blue", "black"))
#par(xpd = TRUE)
lines(x = c(0, 8935), y = c(50, 50), lty = 1, col = "red", lwd = 2)
text(250, 100, "y=50", col = "red")
#abline(h = 50, col = "red")
dev.off()

#remove genes with <50 cells expressing gene, save resulting df.
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180609_CleanUpLawfCluterExpdGenes/")
k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted_top5304 <- k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted[-c(cutoff+1:8935), ]
write.table(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted_top5304, 
            "k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted_top5304.txt", sep = "\t")
top5304genes <- as.vector(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted_top5304$FBN_Flybase)
write(top5304genes, "top5304genes.txt")

#Got GO terms for top5304genes.txt FBgn IDs from flybase: http://flybase.org/batchdownload. Choose "Precomputed files" from data source!
a<-read.csv("top5304genes_goterms.txt", sep = "\t", header = FALSE)
top5304_GOTerms <- data.frame(a$V2, a$V5)
colnames(top5304_GOTerms) <-c("FlyBaseID", "GO_Term")
head(top5304_GOTerms)
unique_GOTERMS<-unique(top5304_GOTerms$GO_Term)
write.table(unique_GOTERMS, "unique_GOTERMS.txt")

# load the GO library, get list of GO ID decriptions
source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
library(GO.db)
goterms <- Term(GOTERM)
goterms <- as.data.frame(goterms)
goterms$GO_ID <- rownames(goterms)
rownames(goterms) <- c()

selectedRows <- (goterms$GO_ID %in% unique_GOTERMS)
goterms_subset <-goterms[selectedRows,]
write.table(goterms_subset, "represented_goterms.txt", sep = "\t")

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180609_CleanUpLawfCluterExpdGenes")
GO_Terms_for_R<-read.csv("GO_Terms_for_R.csv")

#Reduce top5304 to selected GO terms only
M_goterms <- (top5304_GOTerms$GO_Term %in% as.vector(GO_Terms_for_R[,1]))
top5304_GOTerms_M <-top5304_GOTerms[M_goterms,]
AG_goterms <- (top5304_GOTerms$GO_Term %in% as.vector(GO_Terms_for_R[,2]))
top5304_GOTerms_AG <-top5304_GOTerms[AG_goterms,]
CA_goterms <- (top5304_GOTerms$GO_Term %in% as.vector(GO_Terms_for_R[,3]))
top5304_GOTerms_CA <-top5304_GOTerms[CA_goterms,]
#remove duplicate FBgn IDs
M_genes_in_Lawfs <- unique(as.vector(top5304_GOTerms_M$FlyBaseID))
AG_genes_in_Lawfs <- unique(as.vector(top5304_GOTerms_AG$FlyBaseID))
CA_genes_in_Lawfs <- unique(as.vector(top5304_GOTerms_CA$FlyBaseID))
#save files
write(M_genes_in_Lawfs, "M_genes_in_Lawfs.txt")
write(AG_genes_in_Lawfs, "AG_genes_in_Lawfs.txt")
write(CA_genes_in_Lawfs, "CA_genes_in_Lawfs.txt")
