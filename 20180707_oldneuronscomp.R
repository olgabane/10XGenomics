#20180706_oldneuronscomp.R
#Compare oldest neurons in trajectory
#July 6, 2019; July 16, 2019

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

#Load normalized and log transformed gbm. Confirmed is 20180531_GBMpulled.R that this code loads correct normalized, log transformed GBM.
#Some code from 20180531_GBMpulled.R
#loads matrix.mtx, with gene names (rows) and cell barcodes (columns) included. 
library(cellrangerRkit)
gbm <- load_cellranger_matrix("/Users/Olga/Downloads/cellranger/count-Neuron")
#normalize matrix gbm using cell ranger. This also removes genes not exp'd in any cells. 
#convert to full matrix (with 0's).
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
GBM_all <- exprs(gbm_log)
GBM_all <- as.matrix(GBM_all)

#Get cell barcodes for old Lawf neurons. Get GBM with old Lawf neurons only. 
PseudostateDF_OldNeurons_Barcodes <- PseudostateDF_OldNeurons$barcode
old_neuron_positions<-match(PseudostateDF_OldNeurons_Barcodes, colnames(GBM_all))
GBM_oldneurons <- GBM_all[,old_neuron_positions]
#Convert matrix to dataframe (will add colummns)
GBM_oldneurons <-as.data.frame(GBM_oldneurons)

#Delete genes that are not expressed in any of the 496 cells
#sum across rows
sum<-c()
for (i in 1:dim(GBM_oldneurons)[1])
  sum[i] <- sum(GBM_oldneurons[i,])
#add sums in column to dataframe
GBM_oldneurons$RowSum <- sum
length(which(GBM_oldneurons$RowSum == 0)) #3707 genes eliminated (not expressed in any old neurons)
#delete all rows for which RowSum == 0. 
zerorows<-which(GBM_oldneurons$RowSum == 0)
GBM_oldneurons_ExpdGenesOnly <- GBM_oldneurons[-zerorows,]
dim(GBM_oldneurons_ExpdGenesOnly)
#Save file (re-load in R, not excel to preserve all data)
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180706_oldneuronscompare/")
write.table(GBM_oldneurons_ExpdGenesOnly, "GBM_oldneurons_ExpdGenesOnly.txt", sep = "\t")

#Subset dataframe for Lim1+, Lim1- old neurons
Lim1_row <- which(rownames(GBM_oldneurons_ExpdGenesOnly) == "FBgn0026411")
Lim1_neg_col<-which(GBM_oldneurons_ExpdGenesOnly[Lim1_row,] == 0)
Lim1_pos_col<-which(GBM_oldneurons_ExpdGenesOnly[Lim1_row,] != 0)
Lim1_pos_col<-Lim1_pos_col[1:292]
GBM_oldneurons_ExpdGenesOnly_Lim1_neg <- GBM_oldneurons_ExpdGenesOnly[,Lim1_neg_col]
GBM_oldneurons_ExpdGenesOnly_Lim1_pos <- GBM_oldneurons_ExpdGenesOnly[,Lim1_pos_col]

#Add column with expression per cell across row
sum<-c()
for (i in 1:dim(GBM_oldneurons_ExpdGenesOnly_Lim1_neg)[1])
  sum[i] <- sum(GBM_oldneurons_ExpdGenesOnly_Lim1_neg[i,])
GBM_oldneurons_ExpdGenesOnly_Lim1_neg$AvgExprsPerCell <- sum/(dim(GBM_oldneurons_ExpdGenesOnly_Lim1_neg)[2])

sum<-c()
for (i in 1:dim(GBM_oldneurons_ExpdGenesOnly_Lim1_pos)[1])
  sum[i] <- sum(GBM_oldneurons_ExpdGenesOnly_Lim1_pos[i,])
GBM_oldneurons_ExpdGenesOnly_Lim1_pos$AvgExprsPerCell <- sum/(dim(GBM_oldneurons_ExpdGenesOnly_Lim1_pos)[2])

###Identify genes differentially expressed among the two groups of old neurons
relative_expression <- data.frame(GBM_oldneurons_ExpdGenesOnly_Lim1_pos$AvgExprsPerCell, GBM_oldneurons_ExpdGenesOnly_Lim1_neg$AvgExprsPerCell)
rownames(relative_expression) <- rownames(GBM_oldneurons_ExpdGenesOnly_Lim1_pos)
colnames(relative_expression) <- c("AvgExprsPerCell_Lim1_pos", "AvgExprsPerCell_Lim1_neg")
relative_expression$pos_div_neg <-relative_expression[,1]/relative_expression[,2]
#Sort by relative expression
relative_expression_sorted <- relative_expression[order(relative_expression$pos_div_neg),]
#Add columnn characterizing relative expression
LessThan5xDIFF <- which(relative_expression_sorted$pos_div_neg > 0.2 & relative_expression_sorted$pos_div_neg < 5)
Zero <- which(relative_expression_sorted$pos_div_neg == 0)
Infinity <- which(relative_expression_sorted$pos_div_neg == Inf)
Pos_5X <- which(relative_expression_sorted$pos_div_neg >= 5 & relative_expression_sorted$pos_div_neg < Inf)
Neg_5X <- which(relative_expression_sorted$pos_div_neg > 0 & relative_expression_sorted$pos_div_neg <= 0.2)
relative_expression_sorted$Description <- NA
relative_expression_sorted[LessThan5xDIFF,4] <- "Less than 5x diff"
relative_expression_sorted[Zero,4] <- "Zero in Lim1+"
relative_expression_sorted[Infinity,4] <- "Zero in Lim1-"
relative_expression_sorted[Pos_5X,4] <- "5X in Lim1+"
relative_expression_sorted[Neg_5X,4] <- "5X in Lim1-"

#Save resulting table at .txt and .csv
write.csv(relative_expression_sorted, "relative_expression_sorted.csv")
write.table(relative_expression_sorted, "relative_expression_sorted.txt")


