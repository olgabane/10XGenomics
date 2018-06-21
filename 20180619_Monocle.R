#20180619_Monocle.R
#Re-cluster subset using even more PC's from PCA. 
#See if this helps get biologically relevant pseudotrajectory

library(monocle)

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

#Make trajectory using 3-6 PCs
#cluster
for(i in 3:6)
{
    gbm_cds_subset <- reduceDimension(gbm_cds_subset,
                                      max_components = 2,
                                      norm_method = 'log',
                                      num_dim = i,   
                                      reduction_method = 'tSNE',
                                      verbose = T)
    gbm_cds_subset <- clusterCells(gbm_cds_subset, verbose = F)
    
    #Create pseudotrajectory using unsupervised method ('dpFeature') - no prior knowledge of early/late genes.
    #get DE genes
    gbm_cds_subset_expressed_genes <- row.names(subset(fData(gbm_cds_subset),
                                                       num_cells_expressed >= 10)) ####THIS IS A PARAMETER TO CONSIDER CHANGING
    clustering_DEG_genes <-
      differentialGeneTest(gbm_cds_subset[gbm_cds_subset_expressed_genes,],
                           fullModelFormulaStr = '~Cluster',
                           cores = 1)
    
    #order genes
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
    
    #Save CellDataSet
    setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
    Name = paste("gbm_cds_subset_pseudotraj_PC", i, ".RData", sep="")
    saveRDS(gbm_cds_subset, file = Name)
}

#Plots trajectories for PCs 3-6
#Functions tobe used in loop below
GeneByCluster<-function(fname, ...){
  png(file=fname)
  arguments<-list(...)
  Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                        gene_short_name %in% c(arguments)))
  print(plot_genes_jitter(gbm_cds_subset[Genes_of_interest,],
                          grouping = "Cluster",
                          min_expr = 0.1))
  dev.off()
}

GeneByState<-function(fname, ...){
  png(file=fname)
  arguments<-list(...)
  Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                        gene_short_name %in% c(arguments)))
  print(plot_genes_jitter(gbm_cds_subset[Genes_of_interest,],
                    grouping = "State",
                    min_expr = 0.1))
  dev.off()
}

Name <-c("gbm_cds_subset_pseudotraj_PC3.RData", "gbm_cds_subset_pseudotraj_PC4.RData", "gbm_cds_subset_pseudotraj_PC5.RData", "gbm_cds_subset_pseudotraj_PC6.RData")

for(i in 3:6){
    #load CellDatSet
    setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
    gbm_cds_subset<-readRDS(file = Name[i-2])
    #change directory
    Directory = paste("PC", i, sep = "")
    setwd(Directory)
    #plot clusters
    png(file=paste0("PC_", i, ".png"))
    print(plot_cell_clusters(gbm_cds_subset, color_by = 'as.factor(Cluster)'))
    dev.off()
    png(file=paste0("PC_", i, "_markers.png"))
    print(plot_cell_clusters(gbm_cds_subset, markers = c( "Lim1", "hbn", "hth", "eya", "elav", "repo", "ase")))
    dev.off()
    #plot rho, delta
    png(file=paste0("rhodelta_PC_", i, ".png"))
    print(plot_rho_delta(gbm_cds_subset, rho_threshold = 2, delta_threshold = 4)) 
    dev.off()
    #plot trajectories along pseudoaxis by Cluster, Pseudotime, State
    png(file=paste0("Traj_by_cluster_PC_", i, ".png"))
    print(plot_cell_trajectory(gbm_cds_subset, color_by = "Cluster"))
    dev.off()
    png(file=paste0("Traj_by_pseudotime_PC_", i, ".png"))
    print(plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime"))
    dev.off()
    png(file=paste0("Traj_by_state_PC_", i, ".png"))
    print(plot_cell_trajectory(gbm_cds_subset, color_by = "State"))
    dev.off()
    #plot trajectories with markers of interest.
    png(file=paste0("Traj_by_markers_lim1eyahbnhth_PC_", i, ".png"))
    print(plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("Lim1", "eya", "hbn", "hth")))
    dev.off()
    png(file=paste0("Traj_by_markers_elavnSybbrp_PC_", i, ".png"))
    print(plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("elav", "nSyb", "brp")))
    dev.off()
    png(file=paste0("Traj_by_markers_elavrepo_PC_", i, ".png"))
    print(plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("elav", "repo")))
    dev.off()
    png(file=paste0("Traj_by_markers_ase_PC_", i, ".png"))
    print(plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("ase")))
    dev.off()
    #plot genes of interest by cluster, state
    Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                          gene_short_name %in% c("Lim1", "hbn", "hth", "eya")))
    GeneByCluster(paste0("GenebyCluster_LHHE_", i, ".png"), "Lim1", "hbn", "hth", "eya")
    GeneByState(paste0("GenebyState_LHHE_", i, ".png"))
    
    Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                          gene_short_name %in% c("elav", "nSyb", "brp")))
    GeneByCluster(paste0("GenebyCluster_ENB_", i, ".png"))
    GeneByState(paste0("GenebyState_ENB_", i, ".png"))
    
    Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                          gene_short_name %in% c("elav", "repo")))
    GeneByCluster(paste0("GenebyCluster_ER_", i, ".png"))
    GeneByState(paste0("GenebyState_ER_", i, ".png"))
    
    Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                          gene_short_name %in% c("ase")))
    GeneByCluster(paste0("GenebyCluster_A_", i, ".png"))
    GeneByState(paste0("GenebyState_A_", i, ".png"))
}






#encapsulated repetitive __ in funtions to streamline code




