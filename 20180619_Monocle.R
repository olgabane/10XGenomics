#20180619_Monocle.R
#Re-cluster subset using even more PC's from PCA. 
#See if this helps get biologically relevant pseudotrajectory

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180618_Monocle/")
#Load CellDataSet (subsetted and not) created in 20180613_Monocle_largercluster.R
gbm_cds_subset<-readRDS(file = "gbm_cds_subset_1.RData")
gbm_cds<-readRDS(file = "gbm_cds_1.RData")

#Cluster subset using 1-40 PC's
gbm_cds_subset <- detectGenes(gbm_cds_subset, min_expr = 0.1)
fData(gbm_cds_subset)$use_for_ordering <-
  fData(gbm_cds_subset)$num_cells_expressed > 0.05 * ncol(gbm_cds_subset)
quartz("title", 6, 6)
plot_pc_variance_explained(gbm_cds_subset, return_all = F)

#Vary number of PCs used to see how this impacts clustering
for(i in 1:40)
{
gbm_cds_subset <- reduceDimension(gbm_cds_subset,
                                  max_components = 2,
                                  norm_method = 'log',
                                  num_dim = i,   
                                  reduction_method = 'tSNE',
                                  verbose = T)
gbm_cds_subset <- clusterCells(gbm_cds_subset, verbose = F)
png(file=paste0("PC_", i, ".png"))
print(plot_cell_clusters(gbm_cds_subset, color_by = 'as.factor(Cluster)'))
dev.off()
png(file=paste0("PC_", i, "_markers.png"))
print(plot_cell_clusters(gbm_cds_subset, markers = c( "Lim1", "hbn", "hth", "eya", "elav", "repo", "ase")))
dev.off()
}
