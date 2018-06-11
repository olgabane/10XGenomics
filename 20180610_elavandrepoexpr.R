#20180610_elavandrepo
#plots for lab meeting


#elav expression in L's
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180531_GBMpulled/")
k20_Lawf_cluster_GBM_wGeneNames<-read.table("k20_Lawf_cluster_GBM_wGeneNames.txt")

elavrow <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "elav")
quartz("t", 4, 2)
par(mar = c(4,4,2,2))
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[9974, 4:855]), 
     pch=16, cex = 0.5, main = "elav in Lawf cells", xlab = "852 'Lawf' cells", ylab = "Log10(Normalized UMI)", 
     cex.axis = 0.7, cex.lab = 0.7, cex.main = 0.7)

#elav and repo expression in all cells
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180531_GBMpulled/")
