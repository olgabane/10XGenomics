#20180613_Monocle_largercluster.R
#Monocle - arrange cells in pseudotime
#Use progenitor, neuron, repo+ clusters this time
#June 13, 2018, continued June 18

#Load monocle
library(monocle)

#load gene-barcode matrix from cellRanger
library(cellrangerRkit)
cellranger_pipestance_path <- "/Users/Olga/Downloads/cellranger/count-Neuron"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)

#monocle with WHOLE data set
fd <- fData(gbm)
colnames(fd)[2] <- "gene_short_name"
gbm_cds <- newCellDataSet(exprs(gbm),
                          phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                          featureData = new("AnnotatedDataFrame", data = fd),
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
gbm_cds <- estimateSizeFactors(gbm_cds)
gbm_cds <- estimateDispersions(gbm_cds)

#filter
gbm_cds <- detectGenes(gbm_cds, min_expr = 0.1)
print(head(fData(gbm_cds)))
expressed_genes <- row.names(subset(fData(gbm_cds),
                                    num_cells_expressed >= 10))

#cluster cells without gene markers 
disp_table <- dispersionTable(gbm_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
gbm_cds <- setOrderingFilter(gbm_cds, unsup_clustering_genes$gene_id)

#two "intermediate" plots
#genes used for clustering based on expression and dispersion.
quartz("title", 10, 10)
plot_ordering_genes(gbm_cds) #table showing genes used to cluster
#Variance explained by each component
quartz("title", 4, 4)
plot_pc_variance_explained(gbm_cds, return_all = F)

#Plot chosen number of clusters. (Chosen by trial and error). 
#Seems that at this point, Monocle clusters based on tsne space, not PCA, which is unfortunate. 
gbm_cds <- reduceDimension(gbm_cds, max_components = 2, num_dim = 6,  ###could increase # dimensions here based on
                           reduction_method = 'tSNE', verbose = T)
gbm_cds <- clusterCells(gbm_cds, num_clusters = 10)
quartz("title", 4, 4)
plot_cell_clusters(gbm_cds, 1, 2, color = "Cluster", cell_size = 0.5)

#get plot_cell_clusters arguments. Found that I can change point size easily. 
args(plot_cell_clusters)

#plot with markers
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds, 1, 2, markers = c("elav", "repo"), cell_size = 0.5)
#plot with markers
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds, 1, 2, markers = c( "Lim1", "hbn", "hth", "eya"), cell_size = 0.5)
#plot with markers
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds, 1, 2, markers = c( "ase", "dpn", "eya"), cell_size = 0.5)
#plot pros separately because pros range higher; washes out differences in above 3 markers
#plot with markers
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds, 1, 2, markers = c( "pros"), cell_size = 0.5)

#subset CellDataSet to include clusters of interest only
pData(gbm_cds)$RelevantClusters <- (pData(gbm_cds)$Cluster == 1 | pData(gbm_cds)$Cluster == 4 | 
                                       pData(gbm_cds)$Cluster == 7 | pData(gbm_cds)$Cluster == 8)
#verify that this picks out the correct cluster
quartz("title", 4, 4)
plot_cell_clusters(gbm_cds, color_by = 'RelevantClusters')
gbm_cds_subset <- gbm_cds[,pData(gbm_cds)$RelevantClusters]
#verify that all values are now 1, 4, 7, or 8
gbm_cds_subset$Cluster
#2473 cells total
length(gbm_cds_subset$Cluster) 

#Create pseudotrajectory using unsupervised method ('dpFeature') - no prior knowledge of early/late genes.
gbm_cds_subset <- detectGenes(gbm_cds_subset, min_expr = 0.1)
dim(pData(gbm_cds_subset)) 
fData(gbm_cds_subset)$use_for_ordering <-
  fData(gbm_cds_subset)$num_cells_expressed > 0.05 * ncol(gbm_cds_subset)
quartz("title", 6, 6)
plot_pc_variance_explained(gbm_cds_subset, return_all = F)
gbm_cds_subset <- reduceDimension(gbm_cds_subset,
                                  max_components = 2,
                                  norm_method = 'log',
                                  num_dim = 8,   #######THIS number defines how many principal components used.
                                  reduction_method = 'tSNE',
                                  verbose = T)
#Cluster
gbm_cds_subset <- clusterCells(gbm_cds_subset, verbose = F)
#Check that clustering makes sense. 
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds_subset, color_by = 'as.factor(Cluster)')
quartz("markers", 6, 6)
plot_cell_clusters(gbm_cds_subset, markers = c( "Lim1", "hbn", "hth", "eya", "elav", "repo", "ase"))

#plot rho, delta (go back to this)
quartz("title", 6, 6)
plot_rho_delta(gbm_cds_subset, rho_threshold = 2, delta_threshold = 4 ) 

#get DE genes
gbm_cds_subset_expressed_genes <- row.names(subset(fData(gbm_cds_subset),
                                                   num_cells_expressed >= 10)) ####THIS IS A PARAMETER TO CONSIDER CHANGING
clustering_DEG_genes <-
  differentialGeneTest(gbm_cds_subset[gbm_cds_subset_expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 1)

quartz("title", 6, 6)
plot(pData(gbm_cds)$num_genes_expressed)

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

#plot trajectories by Cluster, Pseudotime, State
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Cluster")
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime")
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "State")

#plot trajectories with markers of interest.
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("Lim1", "eya", "hbn", "hth"))
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("elav", "nSyb", "brp"))
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("elav", "repo"))
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("ase"))

#plot gene expression by cluster or state
#I think min_expr = 0.1 represent UMI of 0. 
#They must be using non-logged expression here, and 0.1 was added to al UMI values in the analysis. 
GeneByCluster<-function(x){
  quartz("title", 6, 6)
  plot_genes_jitter(gbm_cds_subset[Genes_of_interest,],
                    grouping = "Cluster",
                    min_expr = 0.1)
}

GeneByState<-function(x){
  quartz("title", 6, 6)
  plot_genes_jitter(gbm_cds_subset[Genes_of_interest,],
                    grouping = "State",
                    min_expr = 0.1)
}

Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                      gene_short_name %in% c("Lim1", "hbn", "hth", "eya")))
GeneByCluster()
GeneByState()

Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                      gene_short_name %in% c("elav", "nSyb", "brp")))
GeneByCluster()
GeneByState()

Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                      gene_short_name %in% c("elav", "repo")))

GeneByCluster()
GeneByState()

Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                      gene_short_name %in% c("ase")))

GeneByCluster()
GeneByState()

#Plot trajectory with gene expression
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("Lim1", "eya", "hbn", "hth"))
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("repo", "elav"))

#Save CellDataSet(subsetted and not). Verified that this is sufficient to re-plot trajectory.
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180618_Monocle/")
saveRDS(gbm_cds_subset, file = "gbm_cds_subset_1.RData")
saveRDS(gbm_cds, file = "gbm_cds_1.RData")
#Load CellDataSet (just as example here, don't run). 
a<-readRDS(file = "gbm_cds_subset_1.RData")




