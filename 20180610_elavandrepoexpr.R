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

#loads matrix.mtx, but with gene names (rows) and cell barcodes (columns) included. Convert to full matrix.
library(cellrangerRkit)
gbm <- load_cellranger_matrix("/Users/Olga/Downloads/cellranger/count-Neuron")
#normalize matrix gbm using cell ranger, name result m2. This also removes genes not exp'd in any cells. 
#convert to full matrix (with 0's).
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
m2 <- exprs(gbm_log)
m2 <- as.matrix(m2)

elav_row <- m2[which(rownames(m2) == "FBgn0260400"), ]
repo_row <- m2[which(rownames(m2) == "FBgn0011701"), ]

elav_repo_df <- data.frame(elav_row, repo_row)
#sort dataframe by elav expression
elav_repo_df_sorted <- elav_repo_df[order(elav_repo_df$elav_row),]

#plot elav vs. repo
quartz("t", 10, 4)
layout(matrix(c(1,2), nrow = 2, byrow = TRUE))
par(mar = c(0,2,0,2), oma = c(4, 1, 1, 1), xpd = NA)
plot(elav_repo_df_sorted$elav_row, pch = 16, col = "red", cex = 0.3, xaxt = "n", xlab = "", cex.axis = 0.7)
plot(elav_repo_df_sorted$repo_row, pch = 16, col = "blue", cex = 0.3, xlab = "7027 cells", cex.axis = 0.7, cex.lab = 0.7)
text(x = -666, y = 1.2, labels = c("Log10(Normalized UMI)"), srt = 90, cex = 0.7)
text(x = -0, y = c(2.4, 1.11), labels = c("elav", "repo"), cex = 0.7, col = c("red", "blue"), font = 2)


do any Lawfs express repo????

