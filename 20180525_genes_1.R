#Use cellrangerRkit to plot gene expression. 
#Will re-do using my G-B Matrix and (made 5.24) and my cluster plots (made 5.22)

#load gene-barcode matrix, normalize, log transform, check dimensions
library(cellrangerRkit)
cellranger_pipestance_path <- "/Users/Olga/Downloads/cellranger/count-Neuron"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180525 Analysis_Gene expression")
genes <-c("elav", "repo", "Lim1", "eya", "hth", "hbn")
tsne_proj <- analysis_results$tsne
quartz("elav-repo-lim1-eya-hth-hbn", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,10))

genes <-c("Pdf", "hbn")
tsne_proj <- analysis_results$tsne


#An attempt to understand limits parameter in visualize_gene_markers() function.
m1_log <- exprs(gbm_log)
m1_log <- as.matrix(m1_log)
range(m1_log)
#result: [1] 0.0000000000000000000 3.1912021200731346404
which(m1_log > 3.19, arr.ind = TRUE) 
#result:              row  col
#result               FBgn0023178 9300 6593
#this is Pdf. So now how 
genes <-c("Pdf")
tsne_proj <- analysis_results$tsne
quartz("pdf1", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,1.5))
quartz("pdf2", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0,10))


