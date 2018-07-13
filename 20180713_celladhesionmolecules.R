#Identify most highly expressed cell adhesion molecules
#July 13, 2018

#Load top expressed genes
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180609_CleanUpLawfCluterExpdGenes/")
top5304 <- read.table("k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted_top5304.txt")

#Load genes with GO term "cell adhesion"
CellAd <- read.table("CA_genes_in_Lawfs.txt")
#clean up rownames
rownames(top5304) <- seq(1:dim(top5304)[1])

#Get rows in top5304 corresponding to FBN IDs in CellAd
rows <- which(match(top5304$FBN_Flybase, CellAd$V1) > 0)
#verify that I get correct number of rows
length(rows) == length(CellAd$V1)

#subset top5304 dataframe for only cell adhesion genes. 
CellAdinLawfs <- top5304[rows, c(1:3,856)]
rownames(CellAdinLawfs) <- seq(1:dim(CellAdinLawfs)[1])

#plot expression of cell adhesion molecules per cell
#Note that values are log10!
quartz("title", 6, 6)
plot(CellAdinLawfs$RowSum/852, pch=16, cex = 0.5)

#Add pseudotrajectory information to expressed cell adhesion molecules
#Get pseudotime DF
setwd("../20180619_Monocle/PC3/AllGenesOverPseudotime/")
PseudotimeDF <- read.csv("AllGenesOverPseudotime_ProcessedJune29.csv")
#Get rows in PseudotimeDF corresponding to FBN IDs in CellAd
#ONLY genes that change over pseudotime trajectory (up or down) are included in PseudotimeDF
rows_P <- which(match(PseudotimeDF$Submitted, CellAd$V1) > 0)
PseudotimeDF_CellAd<-PseudotimeDF[rows_P, ]
rownames(PseudotimeDF_CellAd) <- seq(1:dim(PseudotimeDF_CellAd)[1])
#Get rows in CellAdinLawfs corresponding to cell adhesion molecules that change over time 
rows_PC<-which(match(CellAdinLawfs$FBN_Flybase, PseudotimeDF_CellAd$Submitted_ID) > 0)
#Label all cell adhesion molecules with change in expression over pseudotrajectory
CellAdinLawfs$Trajectory <- NA
a<-match(CellAdinLawfs$FBN_Flybase, PseudotimeDF_CellAd$Submitted_ID)
for (i in 1:length(a))
CellAdinLawfs$Trajectory[i]  = as.character(PseudotimeDF_CellAd$Trajectory[a[i]])
#Label all cell adhesion molecules that don't change as "No_change"
CellAdinLawfs[which(is.na(CellAdinLawfs$Trajectory)), 5] <- "No_Change"

