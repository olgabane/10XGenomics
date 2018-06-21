#20180621_Monocle_4pseudotrajectories.R
#June 21, 2018
#A closer look at 4 pseudotrajectories made using 3-6 principal components

library(monocle)

#Load CellDataSets
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
Name <-c("gbm_cds_subset_pseudotraj_PC3.RData", "gbm_cds_subset_pseudotraj_PC4.RData", "gbm_cds_subset_pseudotraj_PC5.RData", "gbm_cds_subset_pseudotraj_PC6.RData")
for(i in 1:4){assign(paste("gbm_cds_subset_PC_", i+2, sep=""), readRDS(file = Name[i]))}

my_genes <- row.names(subset(fData(gbm_cds_subset_PC_3),
                             gene_short_name %in% c("Lim1", "repo", "elav", "nSyb", "brp", "ase", "eya", "hbn", "hth", "dpr1", "unc-5")))
gbm_cds_subset_PC_3 <- gbm_cds_subset_PC_3[my_genes,]
quartz("title", 6, 6)
plot_genes_in_pseudotime(gbm_cds_subset_PC_3, color_by = "State", nrow = 3, ncol = 4, panel_order = NULL)


my_genes <- row.names(subset(fData(gbm_cds_subset_PC_3),
                             gene_short_name %in% c("nSyb", "brp", "ase")))
gbm_cds_subset_PC_3 <- gbm_cds_subset_PC_3[my_genes,]
quartz("title", 6, 6)
plot_genes_in_pseudotime(gbm_cds_subset_PC_3, color_by = "State")

my_genes <- row.names(subset(fData(gbm_cds_subset_PC_3),
                             gene_short_name %in% c("eya", "hbn", "hth")))
gbm_cds_subset_PC_3 <- gbm_cds_subset_PC_3[my_genes,]
quartz("title", 6, 6)
plot_genes_in_pseudotime(gbm_cds_subset_PC_3, color_by = "State")

#JUNE 22: Group these genes in more logical way
#Plot for PC3-6. Are they similar?
#Choose which to continue with 
