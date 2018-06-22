#20180621_Monocle_4pseudotrajectories.R
#June 21, 2018
#A closer look at 4 pseudotrajectories made using 3-6 principal components

library(monocle)

#Store names of four CellDataSets created using 3-6 principal components in vector
Name <-c("gbm_cds_subset_pseudotraj_PC3.RData", "gbm_cds_subset_pseudotraj_PC4.RData", "gbm_cds_subset_pseudotraj_PC5.RData", "gbm_cds_subset_pseudotraj_PC6.RData")

#Loop over the 4 CellDataSets, plot candidate gene expression as function of pseudotime for pseudotime trajectories made with 3-6 principal components
for(i in 3:6){
#Load CellDataSet
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
gbm_current <- readRDS(file = Name[i-2])
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


