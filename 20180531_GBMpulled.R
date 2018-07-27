#20180531_GBMpulled.R 
#Get GBM that has only cells within clusters pulled yesterday from k=10 and k=20 ("20180530_pull_indiv_clusters.R")
#Some code taken from "20180524_expression_1.R (which has some mistakes in it that have been fixed with this code"

#load sparse matrix (.mtx) file, convert to full matrix (with 0's)
setwd("/Users/Olga/Downloads/cellranger/count-Neuron/outs/filtered_gene_bc_matrices/cellranger/")
library(Matrix)
m1 <- readMM('matrix.mtx')
m1 <- as.matrix(m1)

#loads matrix.mtx, but with gene names (rows) and cell barcodes (columns) included. Don't convert to full matrix (with 0's).
#Will convert after normalization.
library(cellrangerRkit)
gbm <- load_cellranger_matrix("/Users/Olga/Downloads/cellranger/count-Neuron")

#normalize matrix m1 manually 
#sum UMIs for all cells (store in vector x)
x<-c()
for(i in 1:7027)
  x[i]<-sum(m1[,i])
#get median of UMI values
med_x <- median(x)
#divide each value in column by sum of UMIs
for(i in 1:7027)
  for(j in 1:17433)
    m1[j,i] <- m1[j,i]/x[i]
#multiply by median UMI count
m1 <- med_x*(m1)
#delete all rows(genes) which have values of 0
m1 <- m1[which(rowSums(m1) > 0),]
#Log transform normalized matrix (see above for calculation)
m1 <-log10(m1+1)

#normalize matrix gbm using cell ranger, name result m2. This also removes genes not exp'd in any cells. 
#convert to full matrix (with 0's).
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
m2 <- exprs(gbm_bcnorm)
m2 <- exprs(gbm_log)
m2 <- as.matrix(m2)

#verify that the two matrices are equal. #MUST ROUND OFF VALUES; values differ ~15-20 sig digits in.
#note: CANNOT use all() or all.equal() due to different attributes. Must compare values.
a<-which(m1 != m2)
b<-which(m1 == m2)
c<-which(round(m1, digits = 10) != round(m2, digits = 10))
d<-which(round(m1, digits = 10) == round(m2, digits = 10))
sprintf("a = %i, b = %i, c = %i, d = %i, a+b = %i, c+d= %i", length(a), length(b), length(c), length(d), length(a)+length(b), length(c)+length(d))
#second way to verify that matrices are equal.
sum(round(m1, digits = 10) - round(m2, digits = 10))

###NOW GET ONLY GENES FROM CLUSTERS OF INTEREST
#Some of the code below is paired down from "20180530_pull_indiv_clusters.R." All check included in original code.
#load Barcode-Gene-Cluster Color_DataFrame for clusters 2-10 and 11-20
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180522 Analysis_scripts and outputs/20180522_script_2_results/")
df_10 <- read.csv("Barcode-Gene-Cluster-Color_DataFrame_20180522.csv")
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180530_k1120Clusters/Resulting_cluster_plots/")
df_20 <- read.csv("Barcode-Gene-Cluster-Color_DataFrame_20180530.csv")

#Get only relevant columns 
df_10_clust10only <- df_10[, c(2, 3, 4, 23, 24, 34)]
df_20_clust20only <- df_20[, c(2, 3, 4, 23, 24, 34)]

#get only rows corresponding to clusters of interest. Plot to be sure it's correct
#k=10
condition <- df_10_clust10only$kmeans_10_Cluster
k10_yto<-df_10_clust10only[condition == 2 | condition ==3 | condition ==5,]
plot(k10_yto$TSNE.1, k10_yto$TSNE.2, col = paste(k10_yto[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("k=10, pulled clusters", adj = 0.05, line = -2), ylim=c(-50,50), xlim = c(-50, 50))
#k=20
condition <- df_20_clust20only$kmeans_20_Cluster
k20_yto<-df_20_clust20only[condition == 1 | condition ==8 | condition ==6 | condition ==12 |condition ==5,]
plot(k20_yto$TSNE.1, k20_yto$TSNE.2, col = paste(k20_yto[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("k=20, pulled clusters", adj = 0.05, line = -2), ylim=c(-50,50), xlim = c(-50, 50))

#get list of barcodes
k10_barcodes <- as.vector(k10_yto$Barcode)
k20_barcodes <- as.vector(k20_yto$Barcode)

###Extract only genes from relevant clusters from GBM.
#Above, I verified that m1 == m2. Since m2 has barcodes and gene names, I will use m2. 
m3<-m2[,k10_barcodes]
m4 <-m2[,k20_barcodes]
#Verify that the correct number of genes is pulled out.
dim(m3)[2] == length(k10_barcodes) #returns TRUE
dim(m4)[2] == length(k20_barcodes) #returns TRUE

####Now that I've pulled glia, Lawfs, progenitors out of larger matrix, I will split this into 3 matrices for the three diff cell types
#continue working with m4 (K=20 pulled clusters)
#GBM with glia, progenitors and Lawfs
#I will save these files a bit further down the code (when I've added gene names)
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180531_GBMpulled/")
k20_repo_progenitor_lawf_cluster_GBM <-m4
#GBM with glia only
k20_repo_cluster <- k20_yto[which(k20_yto[,6]=="red" | k20_yto[,6]=="darkred"), ]
k20_repo_cluster_barcodes <- as.vector(k20_repo_cluster$Barcode)
k20_repo_cluster_GBM  <-m2[,k20_repo_cluster_barcodes]
#GBM with progenitor only 
k20_progenitor_cluster <- k20_yto[which(k20_yto[,6]=="cyan"), ]
k20_progenitor_cluster_barcodes <- as.vector(k20_progenitor_cluster$Barcode)
k20_progenitor_cluster_GBM  <-m2[,k20_progenitor_cluster_barcodes]
#GBM with Lawfs only
k20_Lawf_cluster <- k20_yto[which(k20_yto[,6]=="pink2" | k20_yto[,6]=="tan3"), ]
k20_Lawf_cluster_barcodes <- as.vector(k20_Lawf_cluster$Barcode)
k20_Lawf_cluster_GBM  <-m2[,k20_Lawf_cluster_barcodes]

#two different ways to verify that split matrices are correct
#check that I have the right number of cells if I add back the three matrices
(dim(k20_Lawf_cluster_GBM)[2] + dim(k20_progenitor_cluster_GBM)[2] 
  + dim(k20_repo_cluster_GBM)[2]) == dim(k20_repo_progenitor_lawf_cluster_GBM)[2] #returns TRUE

#plot all 5 clusters, and repo, progenitor and Lawf clusters broken up
condition <- df_20_clust20only$kmeans_20_Cluster
k20_yto<-df_20_clust20only[condition == 1 | condition ==8 | condition ==6 | condition ==12 |condition ==5,]
quartz("title", 10, 10)
par(mfrow=c(2,2), mar=c(4, 4, 3, 2))
plot(k20_yto$TSNE.1, k20_yto$TSNE.2, col = paste(k20_yto[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, ylim=c(-50,50), xlim = c(-50, 50),
     text(-50, c(45, 40, 35, 30, 25, -40), labels=c("Glia", "Glia", "Progenitors", "Lawfs", "Lawfs", "2510 cells"), 
          col =  c('red', 'darkred', 'cyan', 'pink2', 'tan3', 'black'), adj = -2, pos=4, font = 2))
plot(k20_yto[which(k20_yto[,6]=="red" | k20_yto[,6]=="darkred"), ]$TSNE.1, 
     k20_yto[which(k20_yto[,6]=="red" | k20_yto[,6]=="darkred"), ]$TSNE.2, col = paste(k20_yto[which(k20_yto[,6]=="red" | k20_yto[,6]=="darkred"), ]$kmeans_20_Color), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, ylim=c(-50,50), xlim = c(-50, 50),
     text(-50, c(45, 40, -40), labels=c("Glia", "Glia", "1119 cells"), 
          col =  c('red', 'darkred', 'black'), adj = -2, pos=4, font = 2))
plot(k20_yto[which(k20_yto[,6]=="cyan"), ]$TSNE.1, 
     k20_yto[which(k20_yto[,6]=="cyan"), ]$TSNE.2, col = paste(k20_yto[which(k20_yto[,6]=="cyan"), ]$kmeans_20_Color), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, ylim=c(-50,50), xlim = c(-50, 50),
     text(-50, c(45, -40), labels=c("Progenitors", "539 cells"), 
          col =  c('cyan', 'black'), adj = -2, pos=4, font = 2))
plot(k20_yto[which(k20_yto[,6]=="pink2" | k20_yto[,6]=="tan3"), ]$TSNE.1, 
     k20_yto[which(k20_yto[,6]=="pink2" | k20_yto[,6]=="tan3"), ]$TSNE.2, col = paste(k20_yto[which(k20_yto[,6]=="pink2" | k20_yto[,6]=="tan3"), ]$kmeans_20_Color), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, ylim=c(-50,50), xlim = c(-50, 50), 
     text(-50, c(45, 40, -40), labels=c("Lawf", "Lawf", "852 cells"), 
          col =  c('pink2', 'tan3', 'black'), adj = -2, pos=4, font = 2))
title("k = 20, pulled clusters", outer=TRUE, line = -2)

#individual graphs for powerpoint
quartz("title", 6, 6)
plot(k20_yto$TSNE.1, k20_yto$TSNE.2, col = paste(k20_yto[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, ylim=c(-50,50), xlim = c(-50, 50),
     text(-50, c(45, 40, 35, 30, 25, -40), labels=c("Glia", "Glia", "Progenitors", "Lawfs", "Lawfs", "2510 cells"), 
          col =  c('red', 'darkred', 'cyan', 'pink2', 'tan3', 'black'), adj = -2, pos=4, font = 2))
quartz("title", 6, 6)
plot(k20_yto[which(k20_yto[,6]=="pink2" | k20_yto[,6]=="tan3"), ]$TSNE.1, 
     k20_yto[which(k20_yto[,6]=="pink2" | k20_yto[,6]=="tan3"), ]$TSNE.2, col = paste(k20_yto[which(k20_yto[,6]=="pink2" | k20_yto[,6]=="tan3"), ]$kmeans_20_Color), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, ylim=c(-50,50), xlim = c(-50, 50), 
     text(-50, c(45, 40, -40), labels=c("Lawf", "Lawf", "852 cells"), 
          col =  c('pink2', 'tan3', 'black'), adj = -2, pos=4, font = 2))

#Get FBN ID's, convert to Gene names and symbols
GenesAsFBN_11954<-rownames(k20_Lawf_cluster_GBM)
#converted to gene name and symbols in FlyBase; load resulting file, rename columns
GenesAsFBN_11954_converted_GeneSymbol <- read.csv("GenesAsFBN_11954_converted_GeneSymbol.txt", sep = "\t")
colnames(GenesAsFBN_11954_converted_GeneSymbol) <- c("FBN_Flybase", "GeneName_FlyBase", "GeneSymbol_Flybase")

#Add gene names to GBM's created earlier
library(tibble)
k20_repo_progenitor_lawf_cluster_GBM_wGeneNames <- data.frame(k20_repo_progenitor_lawf_cluster_GBM)
k20_repo_progenitor_lawf_cluster_GBM_wGeneNames <- add_column(k20_repo_progenitor_lawf_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$FBN_Flybase, .before=1)
k20_repo_progenitor_lawf_cluster_GBM_wGeneNames <- add_column(k20_repo_progenitor_lawf_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$GeneName_FlyBase, .before=2)
k20_repo_progenitor_lawf_cluster_GBM_wGeneNames <- add_column(k20_repo_progenitor_lawf_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$GeneSymbol_Flybase, .before=3)
colnames(k20_repo_progenitor_lawf_cluster_GBM_wGeneNames)[1:3] <- c("FBN_Flybase", "GeneName_FlyBase", "GeneSymbol_Flybase")
k20_Lawf_cluster_GBM_wGeneNames <- data.frame(k20_Lawf_cluster_GBM)
k20_Lawf_cluster_GBM_wGeneNames <- add_column(k20_Lawf_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$FBN_Flybase, .before=1)
k20_Lawf_cluster_GBM_wGeneNames <- add_column(k20_Lawf_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$GeneName_FlyBase, .before=2)
k20_Lawf_cluster_GBM_wGeneNames <- add_column(k20_Lawf_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$GeneSymbol_Flybase, .before=3)
colnames(k20_Lawf_cluster_GBM_wGeneNames)[1:3] <- c("FBN_Flybase", "GeneName_FlyBase", "GeneSymbol_Flybase")
k20_repo_cluster_GBM_wGeneNames <- data.frame(k20_repo_cluster_GBM)
k20_repo_cluster_GBM_wGeneNames <- add_column(k20_repo_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$FBN_Flybase, .before=1)
k20_repo_cluster_GBM_wGeneNames <- add_column(k20_repo_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$GeneName_FlyBase, .before=2)
k20_repo_cluster_GBM_wGeneNames <- add_column(k20_repo_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$GeneSymbol_Flybase, .before=3)
colnames(k20_repo_cluster_GBM_wGeneNames)[1:3] <- c("FBN_Flybase", "GeneName_FlyBase", "GeneSymbol_Flybase")
k20_progenitor_cluster_GBM_wGeneNames <- data.frame(k20_progenitor_cluster_GBM)
k20_progenitor_cluster_GBM_wGeneNames <- add_column(k20_progenitor_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$FBN_Flybase, .before=1)
k20_progenitor_cluster_GBM_wGeneNames <- add_column(k20_progenitor_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$GeneName_FlyBase, .before=2)
k20_progenitor_cluster_GBM_wGeneNames <- add_column(k20_progenitor_cluster_GBM_wGeneNames, GenesAsFBN_11954_converted_GeneSymbol$GeneSymbol_Flybase, .before=3)
colnames(k20_progenitor_cluster_GBM_wGeneNames)[1:3] <- c("FBN_Flybase", "GeneName_FlyBase", "GeneSymbol_Flybase")

#Save files
#Note again that .txt files appear to load OK in excel, but some rows are missing. 
#However, I verified in bash using wc that there are 11955 lines, as expected, in each of these files (wc -l .txt)
#SO, I should only load these in R, not in excel!
write.table(k20_repo_cluster_GBM_wGeneNames, "k20_repo_cluster_GBM_wGeneNames.txt", sep = "\t")
write.table(k20_progenitor_cluster_GBM_wGeneNames, "k20_progenitor_cluster_GBM_wGeneNames.txt", sep = "\t")
write.table(k20_Lawf_cluster_GBM_wGeneNames, "k20_Lawf_cluster_GBM_wGeneNames.txt", sep = "\t")
write.table(k20_repo_progenitor_lawf_cluster_GBM_wGeneNames, "k20_repo_progenitor_lawf_cluster_GBM_wGeneNames.txt", sep = '\t')

#Save relevant plots in Demo folder
setwd("/Users/Olga/Desktop/Demo/20180531")
png("PulledClusters.png", 600, 600, res =72)
plot(k20_yto$TSNE.1, k20_yto$TSNE.2, col = paste(k20_yto[,6]), 
     xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.5, title("k=20, pulled clusters", adj = 0.05, line = -2), ylim=c(-50,50), xlim = c(-50, 50))
dev.off()
