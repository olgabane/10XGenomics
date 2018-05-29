#Use cellrangerRkit to plot gene expression. 
#Goal: identify clusters that are elav+ but not Lawf1. 
#Start: May 28, 2018. Continued: May 29.

#load gene-barcode matrix, normalize, log transform, check dimensions
library(cellrangerRkit)
cellranger_pipestance_path <- "/Users/Olga/Downloads/cellranger/count-Neuron"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
analysis_results <- load_cellranger_analysis_results(cellranger_pipestance_path)
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))

#Neuron vs. glia
genes <-c("elav", "repo")
quartz("Neuron vs. glia", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))


#Lawf markers
genes <-c("eya", "hth", "hbn", "Lim1")
quartz("Lawf markers", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Lawf/eg/mg progenitor markers (ase+, pros+, dpn-)
#This makes perfect sense. There are ase+, pros+, dpn- cells inb/w Lawfs and glia, which 
#are the progenitor cells (also marked by 10C12-Gal4)
#Also note that there are dpn+ cells, which cluster near Lawf/egmg progenitors (so these are some NBs)
genes <-c("ase", "pros", "eya", "dpn")
quartz("ase, pros, eya, dpn", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Progenitor should express DE-Cad.
genes <-c("shg")
quartz("shg", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))


#known glial markers
#learned one thing: gcm turns off in Lawfs (reporter stays on, based on Vil's movies)
genes <-c("repo", "Dll", "gcm", "bi")
quartz("Known glial markers", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Xin's table
#all genes in Xin's table (dac is lamina marker, vvl = dfr)
genes <-c("Dll", "repo", "ap", "toy", 
          "D", "oc", "dac", "vvl", "slp1", 
          "run", "ey", "Lim3", "hth", 
          "bsh", "Lim3", "Vsx1", "Vsx2", "svp", "pros")
quartz("All of Xin's genes", 16,10)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))
#only genes significantly expressed in specific clusters 10X data
#NOTE: Dll, ap, ey are in Xavier's paper; non-overlappingly exp'd in medulla neurons
genes <-c("Dll", "repo", "ap", "toy", 
          "oc", "dac", "ey", "hth", 
          "Vsx1", "Vsx2", "pros")
quartz("Xin's genes expd in 10X data", 16,10)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))
#only genes significantly expressed in specific clusters 10X data
#minus Dll, repo (only in glia cluster, as expected)
#re-ordered
genes <-c("toy", "dac", "ey", "pros", "hth",
          "Vsx1", "Vsx2", "ap", "oc")
quartz("Xin's genes expd in unidentified clusters", 16,10)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#####################

#All Drosophila TFs

#Netrin, Ephrin, Semaphorin, Slit receptor
genes <-c("unc-5", "fra", "Eph", "PlexA", "PlexB",
          "robo1", "robo2", "robo3")
quartz("Guidance molecule receptors", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Bsh (Mi1, L4, L5)
genes <-c("bsh")
quartz("bsh", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Spatial factors
genes <-c("Vsx1", "Vsx2", "Optix", "Rx")
quartz("Vsx1, Vsx2, optix, Rx", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#dying cells
genes <-c("rpr", "hid", "grim")
quartz("Rx", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#cadherin
genes <-c("shg")
quartz("shg", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#chinmo 
genes <-c("chinmo")
quartz("chinmo", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#fas2 
genes <-c("Fas2")
quartz("Fas2", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#tOPC signaling molecules
genes <-c("wg", "dpp")
quartz("wg, dpp", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Dscam
genes <-c("Dscam2", "Dscam1")
quartz("Dscam2", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#Dscam2 and Lim1 comparison
genes <-c("Dscam2", "Lim1")
quartz("Dscam2", 8,5)
visualize_gene_markers(gbm_log, genes, tsne_proj[c("TSNE.1", "TSNE.2")], limits = c(0, 1.5))

#### PICK OUT CLUSTER
#### ID genes with least to most variance
###Is lim1 in highest variance??
###If yes, do other genes correlate with Lim1
###Take ALL genes, get GO terms. ID TF's, etc.
###Monocle with progenitor cluster, glia cluster, Lawf cluster

