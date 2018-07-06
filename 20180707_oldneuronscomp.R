#20180706_oldneuronscomp.R
#Compare oldest neurons in trajectory
#July 6, 2019

#Load CellDataSet in which PC=3 clusters were used to make pseudotrajectory
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/")
library(monocle)
gbm_cds_subset_pseudotraj_PC3 = readRDS('gbm_cds_subset_pseudotraj_PC3.RData')
head(pData(gbm_cds_subset_pseudotraj_PC3))

#Get only barcode, Cluster, Pseudotime, State data
PseudostateDF <- data.frame(gbm_cds_subset_pseudotraj_PC3$barcode, gbm_cds_subset_pseudotraj_PC3$Cluster, 
                            gbm_cds_subset_pseudotraj_PC3$State, gbm_cds_subset_pseudotraj_PC3$Pseudotime)
colnames(PseudostateDF)<-c("barcode", "Cluster", "State", "Pseudotime")

#Isolate oldest neurons (Pseudotime >= 30)
PseudostateDF_OldNeurons  <- subset(PseudostateDF, Pseudotime >=30)

