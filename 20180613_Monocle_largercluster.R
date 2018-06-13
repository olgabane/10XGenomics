#20180613_Monocle_largercluster.R
#Monocle - arrange cells in pseudotime
#Use progenitor, neuron, repo+ clusters this time
#June 13, 2018

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
gbm_cds <- reduceDimension(gbm_cds, max_components = 2, num_dim = 6,
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
#First, need to add column to phenotypic data 
###THIS DOESN'T WORK!! Will have to subset more manually. 
pData(gbm_cds)$RelevantClusters <- (pData(gbm_cds)$Cluster == 1 || pData(gbm_cds)$Cluster == 5 || 
  pData(gbm_cds)$Cluster == 7 || pData(gbm_cds)$Cluster == 8)
#verify that this picks out the correct cluster
quartz("title", 4, 4)
plot_cell_clusters(gbm_cds, color_by = 'RelevantClusters')
gbm_cds_subset <- gbm_cds[,pData(gbm_cds)$RelevantClusters]
#verify that all values are now 1
gbm_cds_subset$Cluster
length(gbm_cds_subset$Cluster) 


length(which(gbm_cds$Cluster == 1)) 

quartz("title", 6, 6)
plot(pData(gbm_cds)$num_genes_expressed)







########NEED TO LOOK AT CELLRANGER OUTPUT: 
#Need to see what each PCA is based on. i.e if some based on cell cycle, or cell stress. (See Seurat FAQ)
#-DO this in Seurat instead? Seurat provides function to do this. 