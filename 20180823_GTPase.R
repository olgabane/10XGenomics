#20180823_GTPase.R
#Goal: GTPase expression
#August 23

#load gene-barcode matrix, normalize, log transform, check dimensions
library(cellrangerRkit)
cellranger_pipestance_path <- "/Users/Olga/Downloads/cellranger/count-Neuron"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))
tsne_proj <- analysis_results$tsne

genes <-c("Rac1", "Rac2", "Mtl")
quartz("title", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 4))

genes <-c("Rho1")
quartz("title", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 4))