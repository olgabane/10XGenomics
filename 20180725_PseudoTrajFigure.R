#20180725_PseudoTrajFigure.R
#July 25, 2018

library(monocle)

#Store names of four CellDataSets created using 3-6 principal components in vector
Name <-c("gbm_cds_subset_pseudotraj_PC3.RData", "gbm_cds_subset_pseudotraj_PC4.RData", "gbm_cds_subset_pseudotraj_PC5.RData", "gbm_cds_subset_pseudotraj_PC6.RData")

#Loop over the 4 CellDataSets, plot candidate gene expression as function of pseudotime for pseudotime trajectories made with 3-6 principal components
#Load CellDataSet for PC = 3
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
gbm_current <- readRDS(file = "gbm_cds_subset_pseudotraj_PC3.RData")
setwd("PC3/")
#plot genes over pseudotime
my_genes <- row.names(subset(fData(gbm_current),
                               gene_short_name %in% c("Lim1", "repo", "elav", "elav", 
                                                      "nSyb", "brp", "ase", "eya", "hbn", "hth", "dpr1", "Frq1", "mbl")))
gbm_current <- gbm_current[my_genes,]
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/")
png(file="", width = 1000, height = 600)
print(plot_genes_in_pseudotime(gbm_current, color_by = "State", nrow = 4, ncol = 4, 
                                 panel_order = c( "repo", "elav", "nSyb", "brp", "Lim1", 
                                                  "eya", "hbn", "hth", "ase", "pros", "dpr1", "Frq1", "mbl")))
dev.off()
}