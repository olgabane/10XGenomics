#20180706_oldneuronscomp.R
#Compare oldest neurons in trajectory
#July 6, 2019; July 16, 2019

#Load CellDataSet in which PC=3 clusters were used to make pseudotrajectory
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
library(monocle)
gbm_cds_subset_pseudotraj_PC3 = readRDS('gbm_cds_subset_pseudotraj_PC3.RData')
head(pData(gbm_cds_subset_pseudotraj_PC3))

#Get only barcode, Cluster, Pseudotime, State data
PseudostateDF <- data.frame(gbm_cds_subset_pseudotraj_PC3$barcode, gbm_cds_subset_pseudotraj_PC3$Cluster, 
                            gbm_cds_subset_pseudotraj_PC3$State, gbm_cds_subset_pseudotraj_PC3$Pseudotime)
colnames(PseudostateDF)<-c("barcode", "Cluster", "State", "Pseudotime")

#Isolate oldest neurons (Pseudotime >= 30)
PseudostateDF_OldNeurons  <- subset(PseudostateDF, Pseudotime >=30)

#Load normalized and log transformed gbm. Confirmed is 20180531_GBMpulled.R that this code loads correct normalized, log transformed GBM.
#Some code from 20180531_GBMpulled.R
#loads matrix.mtx, with gene names (rows) and cell barcodes (columns) included. 
library(cellrangerRkit)
gbm <- load_cellranger_matrix("/Users/Olga/Downloads/cellranger/count-Neuron")
#normalize matrix gbm using cell ranger. This also removes genes not exp'd in any cells. 
#convert to full matrix (with 0's).
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
GBM_all <- exprs(gbm_log)
GBM_all <- as.matrix(GBM_all)

#Get cell barcodes for old Lawf neurons. Get GBM with old Lawf neurons only. 
PseudostateDF_OldNeurons_Barcodes <- PseudostateDF_OldNeurons$barcode
old_neuron_positions<-match(PseudostateDF_OldNeurons_Barcodes, colnames(GBM_all))
GBM_oldneurons <- GBM_all[,old_neuron_positions]


