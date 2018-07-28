#20180727: A synthesis of code from the following scripts:
#20180613_Monocle_largercluster.R
#20180619_Monocle.R
#20180621__Monocle_4pseudotrajectories.R
#Goal: Use Monocle to arrange cells clusters pulled in pseudotime trajecotry

#Load monocle
library(monocle)

#load gene-barcode matrix from cellRanger. Convert to CellDataSet for Monocle.
library(cellrangerRkit)
cellranger_pipestance_path <- "/Users/Olga/Downloads/cellranger/count-Neuron"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
fd <- fData(gbm)
colnames(fd)[2] <- "gene_short_name"
gbm_cds <- newCellDataSet(exprs(gbm),
                          phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                          featureData = new("AnnotatedDataFrame", data = fd),
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())

#Calculate two pieces of information about the data (because modeled as negative binomial distribution):
#Neg binomial distribution: discrete event (count belongs to single gene of many). Low probability descrete event
#but does not follow Poisson distibution because variability in read counts too large (in P, variance = mean).
#in sequencing data, variance > mean. 
#SizeFactor is used to normalize for library size in each cell (UMI count/mRNA recovered)
#Becuase of overdispersion in sequencing data, use neg binomial distributaion --> requires dispersion parameter to estimate variance correctly.
gbm_cds <- estimateSizeFactors(gbm_cds)
gbm_cds <- estimateDispersions(gbm_cds)

#Filter low quality data
#Calculate number of cells that express gene, min_expr defined
gbm_cds <- detectGenes(gbm_cds, min_expr = 0.1)
print(head(fData(gbm_cds)))
#vector of genes expressed in at least 10 cells of data set. To be used in pseudotrajectory.
expressed_genes <- row.names(subset(fData(gbm_cds),
                                    num_cells_expressed >= 10))

#cluster cells (unsupervised, no gene markers) 
#choose which cells to use for clustering. Choose genes that vary and are expressed at some threshold.
#setOrdering filter marks genes to be used for clustering. 
disp_table <- dispersionTable(gbm_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
gbm_cds <- setOrderingFilter(gbm_cds, unsup_clustering_genes$gene_id)
#Reduce dimensions using tsne. num_dim is number PCA components used to reduce dimensions. Determined by plotting variance explained by each component.
#clusterCells() uses density peak clustering rather than k-means.
gbm_cds <- reduceDimension(gbm_cds, max_components = 2, num_dim = 6,  
                           reduction_method = 'tSNE', verbose = T)
gbm_cds <- clusterCells(gbm_cds, num_clusters = 10)
quartz("title", 4, 4)
plot_cell_clusters(gbm_cds, 1, 2, color = "Cluster", cell_size = 0.5)

#plot with markers to identify clusters of interest
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds, 1, 2, markers = c("elav", "repo",  "Lim1", "hbn", "hth", "eya", "ase", "dpn"), cell_size = 0.5)

#subset CellDataSet to include clusters of interest only
pData(gbm_cds)$RelevantClusters <- (pData(gbm_cds)$Cluster == 1 | pData(gbm_cds)$Cluster == 4 | 
                                      pData(gbm_cds)$Cluster == 7 | pData(gbm_cds)$Cluster == 8)

#verify that this picks out the correct cluster
quartz("title", 4, 4)
plot_cell_clusters(gbm_cds, color_by = 'RelevantClusters')
gbm_cds_subset <- gbm_cds[,pData(gbm_cds)$RelevantClusters]
#verify that all values are now only clusters pulled out
gbm_cds_subset$Cluster
#2473 cells total (varies because tsne and clustering are stochastic)
length(gbm_cds_subset$Cluster) 

#Make trajectory using 3-6 PCs. ML algorithm: reverse graph embedding.
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
  #Choose which genes to use: get genes differentially expressed between clusters. 
  gbm_cds_subset_expressed_genes <- row.names(subset(fData(gbm_cds_subset),
                                                     num_cells_expressed >= 10)) ####THIS IS A PARAMETER TO CONSIDER CHANGING
  clustering_DEG_genes <-
    differentialGeneTest(gbm_cds_subset[gbm_cds_subset_expressed_genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = 1)
  
  #select top 1000 ordering genes to order cells 
  gbm_cds_subset_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  
  gbm_cds_subset <- setOrderingFilter(gbm_cds_subset, ordering_genes = gbm_cds_subset_ordering_genes)

  gbm_cds_subset<- reduceDimension(gbm_cds_subset, method = 'DDRTree')
  
  gbm_cds_subset <- orderCells(gbm_cds_subset)
  
  #Save CellDataSet
  setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
  Name = paste("gbm_cds_subset_pseudotraj_PC", i, ".RData", sep="")
  saveRDS(gbm_cds_subset, file = Name)
}

#Plots trajectories for PCs 3-6
#Functions to be used in loop below
GeneByCluster<-function(fname, ...){
  png(file=fname)
  arguments<-list(...)
  Genes_of_interest <- row.names(subset(fData(gbm_cds_subset), gene_short_name %in% c(arguments)))
  print(plot_genes_jitter(gbm_cds_subset[Genes_of_interest,], grouping = "Cluster", min_expr = 0.1))
  dev.off()
}

GeneByState<-function(fname, ...){
  png(file=fname)
  arguments<-list(...)
  Genes_of_interest <- row.names(subset(fData(gbm_cds_subset), gene_short_name %in% c(arguments)))
  print(plot_genes_jitter(gbm_cds_subset[Genes_of_interest,], grouping = "State", min_expr = 0.1))
  dev.off()
}

PlotMarkersOverTraj <- function(fname, ...){
  png(file=fname)
  arguments <- list(...)
  print(plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c(arguments)))
  dev.off()
}

PlotTrajectory <- function(fname, how){
  png(file=fname)
  print(plot_cell_trajectory(gbm_cds_subset, color_by = how))
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
  PlotTrajectory(paste0("Traj_by_cluster_PC_", i, ".png"), "Cluster")
  PlotTrajectory(paste0("Traj_by_pseudotime_PC_", i, ".png"), "Pseudotime")
  PlotTrajectory(paste0("Traj_by_state_PC_", i, ".png"), "State")
  
  #plot trajectories with markers of interest.
  PlotMarkersOverTraj(paste0("Traj_by_markers_lim1eyahbnhth_PC_", i, ".png"), "Lim1", "eya", "hbn", "hth")
  PlotMarkersOverTraj(paste0("Traj_by_markers_elavnSybbrp_PC_", i, ".png"), "elav", "nSyb", "brp")
  PlotMarkersOverTraj(paste0("Traj_by_markers_elavrepo_PC_", i, ".png"), "elav", "repo")
  PlotMarkersOverTraj(paste0("Traj_by_markers_ase_PC_", i, ".png"), "ase")
  
  #plot genes of interest by cluster, state
  GeneByCluster(paste0("GenebyCluster_LHHE_", i, ".png"), "Lim1", "hbn", "hth", "eya")
  GeneByState(paste0("GenebyState_LHHE_", i, ".png"), "Lim1", "hbn", "hth", "eya")
  GeneByCluster(paste0("GenebyCluster_ENB_", i, ".png"), "elav", "nSyb", "brp")
  GeneByState(paste0("GenebyState_ENB_", i, ".png"), "elav", "nSyb", "brp")
  GeneByCluster(paste0("GenebyCluster_ER_", i, ".png"), "elav", "repo")
  GeneByState(paste0("GenebyState_ER_", i, ".png"), "elav", "repo")
  GeneByCluster(paste0("GenebyCluster_A_", i, ".png"), "ase")
  GeneByState(paste0("GenebyState_A_", i, ".png"), "ase")
}

#Store names of four CellDataSets created using 3-6 principal components in vector
Name <-c("gbm_cds_subset_pseudotraj_PC3.RData", "gbm_cds_subset_pseudotraj_PC4.RData", "gbm_cds_subset_pseudotraj_PC5.RData", "gbm_cds_subset_pseudotraj_PC6.RData")

#Loop over the 4 CellDataSets, plot candidate gene expression as function of pseudotime for pseudotime trajectories made with 3-6 principal components
for(i in 3:6){
  #Load CellDataSet
  setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
  gbm_current <- readRDS(file = Name[i-2])
  setwd("/Users/Olga/Desktop/Demo/20180727/")
  setwd(paste("PC", i, sep = ""))
  #plot genes over pseudotime
  my_genes <- row.names(subset(fData(gbm_current),
                               gene_short_name %in% c("Lim1", "repo", "elav", "elav", 
                                                      "nSyb", "brp", "ase", "eya", "hbn", "hth", "dpr1", "Frq1", "mbl")))
  gbm_current <- gbm_current[my_genes,]
  png(file=paste0("GenesOverPseudotime_PC_", i, ".png"), width = 1000, height = 600)
  print(plot_genes_in_pseudotime(gbm_current, color_by = "State", nrow = 4, ncol = 4, 
                                 panel_order = c( "repo", "elav", "nSyb", "brp", "Lim1", 
                                                  "eya", "hbn", "hth", "ase", "pros", "dpr1", "Frq1", "mbl")))
  dev.off()
}

#Goal: ID genes with specific, correlated temporal expression patterns
#Choose single CellDataSet (made using 3, 4, 5 or 6 PCs), plot ALL genes expressed in >10 cells as function of pseudotime
#For now, I will use CellDataSet in which trajectory was made using 3 PCs
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
gbm_current <- readRDS(file = "gbm_cds_subset_pseudotraj_PC3.RData")
length(row.names(subset(fData(gbm_current)))) #DF has all 17433 genes, whether or not they are expressed. 
#Get DF with genes expressed in at least 10 cells
gbm_current_expressedin10 <- fData(gbm_current)[which (fData(gbm_current)$num_cells_expressed >= 10),] 
dim(gbm_current_expressedin10) #Leaves 8031 genes

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/PC3/AllGenesOverPseudotime/")
for (i in 1:dim(gbm_current_expressedin10)[1])
{
  my_genes <- row.names(gbm_current_expressedin10)
  my_genes_current <- my_genes[i]
  gbm_current_singlegene <- gbm_current[my_genes_current,]
  png(file=paste0(my_genes_current, ".png"), width = 1000, height = 600)
  print(plot_genes_in_pseudotime(gbm_current_singlegene, color_by = "State"))
  dev.off()
}



