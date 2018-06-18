#20180618_Monocle.R
#Re-cluster subset using more PC's from PCA. 
#See if this helps get biologically relevant pseudotrajectory

library(monocle)

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180618_Monocle/")
#Load CellDataSet (subsetted and not) created in 20180613_Monocle_largercluster.R
gbm_cds_subset<-readRDS(file = "gbm_cds_subset_1.RData")
gbm_cds<-readRDS(file = "gbm_cds_1.RData")

#Cluster subset using first 20 PC's (previously: used 8)
gbm_cds_subset <- detectGenes(gbm_cds_subset, min_expr = 0.1)
fData(gbm_cds_subset)$use_for_ordering <-
  fData(gbm_cds_subset)$num_cells_expressed > 0.05 * ncol(gbm_cds_subset)
quartz("title", 6, 6)
plot_pc_variance_explained(gbm_cds_subset, return_all = F)
gbm_cds_subset <- reduceDimension(gbm_cds_subset,
                                  max_components = 2,
                                  norm_method = 'log',
                                  num_dim = 20,   #######THIS number defines how many principal components used.
                                  reduction_method = 'tSNE',
                                  verbose = T)
gbm_cds_subset <- clusterCells(gbm_cds_subset, verbose = F)
#Check that clustering makes sense. 
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds_subset, color_by = 'as.factor(Cluster)')
quartz("markers", 6, 6)
plot_cell_clusters(gbm_cds_subset, markers = c( "Lim1", "hbn", "hth", "eya", "elav", "repo", "ase"))

#plot rho, delta (go back to this)
quartz("title", 6, 6)
plot_rho_delta(gbm_cds_subset, rho_threshold = 2, delta_threshold = 4 ) 

#Create pseudotrajectory using unsupervised method ('dpFeature') - no prior knowledge of early/late genes.
#get DE genes
gbm_cds_subset_expressed_genes <- row.names(subset(fData(gbm_cds_subset),
                                                   num_cells_expressed >= 10)) ####THIS IS A PARAMETER TO CONSIDER CHANGING
clustering_DEG_genes <-
  differentialGeneTest(gbm_cds_subset[gbm_cds_subset_expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 1)

#plot along pseudoaxis
gbm_cds_subset_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

gbm_cds_subset <-
  setOrderingFilter(gbm_cds_subset,
                    ordering_genes = gbm_cds_subset_ordering_genes)
plot_ordering_genes(gbm_cds_subset)

gbm_cds_subset<-
  reduceDimension(gbm_cds_subset, method = 'DDRTree')

gbm_cds_subset <-
  orderCells(gbm_cds_subset)

