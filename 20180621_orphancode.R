#20180621_orphancode.R
#Unused code pieces that I want to keep for future reference.

#save all open plots
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
for(d in dev.list()) {
  dev.set(d)
  Name = paste("Plot", d, "_PC3.png", sep="")
  dev.copy(png, Name)
  dev.off()
}

#Save CellDataSet
saveRDS(gbm_cds_subset, file = "gbm_cds_subset_PC3.RData")
#Load CellDataSet (just as example here, don't run). 
a<-readRDS(file = "gbm_cds_subset_1.RData")

#this function is worth looking at further. 
HSMM_expressed_genes <-  row.names(subset(fData(HSMM_myo),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Hours")

#get pseudotime axis. Still unclear how to get tsne plot coordinates. 
pData(gbm_cds_subset)$Pseudotime
