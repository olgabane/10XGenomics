#20180608_Monocle.R
#Monocle - arrange cells in pseudotime
#June 8, 2018


#Install Bioconductor, Monocle and verify
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("monocle")
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

gbm_cds <- detectGenes(gbm_cds, min_expr = 0.1)
print(head(fData(gbm_cds)))
expressed_genes <- row.names(subset(fData(gbm_cds),
                                    num_cells_expressed >= 10))

#cluster cells without gene markers 
disp_table <- dispersionTable(gbm_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
gbm_cds <- setOrderingFilter(gbm_cds, unsup_clustering_genes$gene_id)
quartz("title", 10, 10)
plot_ordering_genes(gbm_cds) #table showing genes used to cluster

quartz("title", 4, 4)
plot_pc_variance_explained(gbm_cds, return_all = F) 
#Plot 10 clusters. Note: I chose 10 by trial and error (captures L cluster well)
gbm_cds <- reduceDimension(gbm_cds, max_components = 2, num_dim = 6,
                           reduction_method = 'tSNE', verbose = T)
gbm_cds <- clusterCells(gbm_cds, num_clusters = 10)
quartz("title", 4, 4)
plot_cell_clusters(gbm_cds, 1, 2, color = "Cluster")

#plot with markers
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds, 1, 2, color = "Cluster", markers = c("elav", "repo"))
#plot with markers
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds, 1, 2, color = "Cluster", markers = c( "Lim1", "hbn", "hth", "eya"))

length(which(gbm_cds$Cluster == 1)) 

#subset CellDataSet to include cluster 1 only; Will use this subset to create pseudotime axis
#First, need to add column to phenotypic data 
pData(gbm_cds)$LawfCluster <- pData(gbm_cds)$Cluster == 1
#verify that this picks out the correct cluster
quartz("title", 4, 4)
plot_cell_clusters(gbm_cds, color_by = 'LawfCluster')
gbm_cds_subset <- gbm_cds[,pData(gbm_cds)$LawfCluster]
#verify that all values are now 1
gbm_cds_subset$Cluster
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
                                  num_dim = 3,
                                  reduction_method = 'tSNE',
                                  verbose = T)
gbm_cds_subset <- clusterCells(gbm_cds_subset, verbose = F)
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds_subset, color_by = 'as.factor(Cluster)')
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds_subset, markers = c( "Lim1", "hbn", "hth", "eya"))
quartz("title", 6, 6)
plot_cell_clusters(gbm_cds_subset, markers = c( "elav"))
#"Can also provide the decision plot for users to check the Ρ, Δ for each cell and decide the threshold for defining the cell clusters".
#keep this in mind, but it's not something I did this first time
quartz("title", 6, 6)
plot_rho_delta(gbm_cds_subset, rho_threshold = 2, delta_threshold = 4 ) ###may need to change these parameters

#get DE genes
gbm_cds_subset_expressed_genes <- row.names(subset(fData(gbm_cds_subset),
                                                   num_cells_expressed >= 10))
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


#plot trajectories in several ways
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Cluster")

quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime")

quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "State")

quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("Lim1", "eya", "hbn", "hth"))
quartz("title", 6, 6)
plot_cell_trajectory(gbm_cds_subset, color_by = "Pseudotime", markers = c("elav", "nSyb", "brp"))

State_Barcode_df <- data.frame(pData(gbm_cds_subset)$barcode, pData(gbm_cds_subset)$State)

Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                      gene_short_name %in% c("Lim1", "hbn", "hth", "eya")))
quartz("title", 6, 6)
plot_genes_jitter(gbm_cds_subset[Genes_of_interest,],
                  grouping = "Cluster",
                  min_expr = 0.1)

Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                      gene_short_name %in% c("Lim1", "hbn", "hth", "eya")))
quartz("title", 6, 6)
plot_genes_jitter(gbm_cds_subset[Genes_of_interest,],
                  grouping = "State",
                  min_expr = 0.1)

Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                      gene_short_name %in% c("elav", "nSyb", "brp")))
quartz("title", 6, 6)
plot_genes_jitter(gbm_cds_subset[Genes_of_interest,],
                  grouping = "Cluster",
                  min_expr = 0.1)

Genes_of_interest <- row.names(subset(fData(gbm_cds_subset),
                                      gene_short_name %in% c("elav", "nSyb", "brp")))
quartz("title", 6, 6)
plot_genes_jitter(gbm_cds_subset[Genes_of_interest,],
                  grouping = "State",
                  min_expr = 0.1)

quartz("title", 6, 6)
plot_cell_clusters(gbm_cds, 1, 2, color = "Cluster", markers = "Lim1")
