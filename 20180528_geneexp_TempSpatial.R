#Use cellrangerRkit to plot gene expression. 
#Goal: identify clusters that are elav+ but not Lawf1. 

#load gene-barcode matrix, normalize, log transform, check dimensions
library(cellrangerRkit)
cellranger_pipestance_path <- "/Users/Olga/Downloads/cellranger/count-Neuron"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))

#ey
genes <-c("ey")
tsne_proj <- analysis_results$tsne
quartz("ey", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Ap
genes <-c("ap")
quartz("ap", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Slp
genes <-c("slp1", "slp2")
quartz("Slp1, slp2", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#D
genes <-c("D")
quartz("D", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Lim3
genes <-c("Lim3")
quartz("Lim3", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Bsh
genes <-c("bsh")
quartz("bsh", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Toy
genes <-c("toy")
quartz("toy", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Dfr
genes <-c("vvl")
quartz("vvl/Dfr", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Hth
genes <-c("hth")
quartz("hth", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Dll
genes <-c("Dll")
quartz("Dll", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#vsx
genes <-c("Vsx1", "Vsx2")
quartz("Vsx1, Vsx2", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#optix
genes <-c("Optix")
quartz("Optix", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#rx
genes <-c("Rx")
quartz("Rx", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#L1-L5:

#dying cells
