#20180601_LawfClusterCloserLook.R
#June 1, 2018

#Load GBM file created 20180531 (see 20180531_GBMpulled.R)
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180531_GBMpulled/")
k20_Lawf_cluster_GBM_wGeneNames<-read.table("k20_Lawf_cluster_GBM_wGeneNames.txt")

#Plot Normalized UMI values 
gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "Lim1")
quartz("title", 10, 10)
layout(matrix(c(1:8), nrow = 4, byrow = TRUE), widths=c(4,1))
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main="Lim1", xlab = "852 'Lawf'cells", ylab = "Log10(Normalized UMI)")
p2=barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0))), main = "Lim1", col = c("red", "blue"))
axis(1, at = p2, labels = c("0", ">0"))

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "hth")
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main="hth", xlab = "852 'Lawf'cells", ylab = "Log10(Normalized UMI)")
p2=barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0))), main = "hth", col = c("red", "blue"))
axis(1, at = p2, labels = c("0", ">0"))

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "hbn")
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main="hbn", xlab = "852 'Lawf'cells", ylab = "Log10(Normalized UMI)")
p2=barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0))), main = "hbn", col = c("red", "blue"))
axis(1, at = p2, labels = c("0", ">0"))

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "eya")
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main="eya", xlab = "852 'Lawf'cells", ylab = "Log10(Normalized UMI)")
p2=barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0))), main = "eya", col = c("red", "blue"))
axis(1, at = p2, labels = c("0", ">0"))



####NEXT: Delete genes with 0 when summed across all cells
  #how many are left?
#####Sum UMIs of all genes. 
  #Look at most highly exp'd
#####Look at most variable genes across Lawfs. 
  #Is Lim1 one of them?
  #If yes, do any correlate w Lim1?



##############

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "elav")
quartz("title", 10, 10)
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main)

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "repo")
quartz("title", 10, 10)
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main)
a<-which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)
b<-which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0)
length(a)
length(b)
#### delete these if they don't have elav in them!!

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "ase")
quartz("title", 10, 10)
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main = "Lim1 in Lawf cluster")
a<-which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)
b<-which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0)