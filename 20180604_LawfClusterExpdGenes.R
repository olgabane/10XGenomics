#20180604_LawfClusterExpdGenes.R
#June 4, 2018

#Load GBM file created 20180531 (see 20180531_GBMpulled.R)
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180531_GBMpulled/")
k20_Lawf_cluster_GBM_wGeneNames<-read.table("k20_Lawf_cluster_GBM_wGeneNames.txt")

###Delete genes that are not expressed in any of the 852 cells
#sum across rows
sum<-c()
for (i in 1:dim(k20_Lawf_cluster_GBM_wGeneNames)[1])
  sum[i] <- sum(k20_Lawf_cluster_GBM_wGeneNames[i, 4:855])
#add sums in column to dataframe
k20_Lawf_cluster_GBM_wGeneNames$RowSum <- sum
length(which(k20_Lawf_cluster_GBM_wGeneNames$RowSum == 0)) #3019 genes eliminated (not expressed in any cells in this cluster)
#delete all rows for which RowSum == 0. I did some manual check to make sure this worked properly.
zerorows<-which(k20_Lawf_cluster_GBM_wGeneNames$RowSum == 0)
k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly <- k20_Lawf_cluster_GBM_wGeneNames[-zerorows,]
dim(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly)
#Save file
#Note again that .txt files appear to load OK in excel, but some rows are missing (see 20180531_GBMpulled.R for details)
#SO, I should only load these in R, not in excel!
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180601_LawfClustersCloserLook/")
write.table(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly, "k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly.txt", sep = "\t")
#Sort by expression (RowSums), in descending order
k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted <-
  k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[order(-k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly$RowSum),]
write.table(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted, "k20_repo_cluster_GBM_wGeneNames_ExpdGenesOnly_sorted.txt", sep = "\t")

###Identify most variable genes
#calculate standard deviation across rows, for expressed genes only
#####..........Neeed to calculate COEFFICIENT OF VARIATION here, b/c higher value will result in higher SD (not normalized)
#####..........COEFFICIENT OF VARIATION = std/mean
#sort based on coefficient of varition
stddev<-c()
for (i in 1:dim(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly)[1])
  stddev[i] <- sd(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[i, 4:855])
mean <-c()
for (i in 1:dim(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly)[1])
  mean[i] <- mean(as.numeric(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[i, 4:855]))
coef_of_var = stddev/mean

k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_SortedbyCOV <- k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly
k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_SortedbyCOV$COV <- coef_of_var
k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_SortedbyCOV <- 
  k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_SortedbyCOV[order(-k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_SortedbyCOV$COV),]
write.table(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_SortedbyCOV, "k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly_SortedbyCOV.txt", sep = "\t")
###THIS DIDN'T WORK TOO WELL...super lowly exp'd genes or ones exp'd in just a few cells have tiny mean, huge correlation value.
  ###Go back and look at STDEV. Lim1 comes out on top there...MBL looks very promising!!
###OR: do this. sort by column: Lim1+ vs Lim1-. Then, sum expr off all other genes in Lim1+ minus Lim1-. Look for biggest difference,.

###Sort by Lim1. PLot ALL expressed genes.
#Not very efficient to keep re-sorting (rather than sort whole df by Lim1 first), but less likely to make mistakes this way.
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180601_LawfClustersCloserLook/Rplots_8935Genes/")
Lim1_row <- which(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly$GeneSymbol_Flybase == "Lim1")
for (i in 1:dim(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly)[1]) 
{
  twogenes_df <- data.frame(as.numeric(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lim1_row, 4:855]), 
                            as.numeric(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[i, 4:855]))
  colnames(twogenes_df) <- c("Lim1", "gene2")
  #reorder based on Lim1 expression 
  twogenes_df_ordered<-twogenes_df[order(twogenes_df$Lim1),]
  #plot
  png(paste(rownames(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[i,]), ".png", sep = ""))
  layout(matrix(c(1,2), nrow = 2, byrow = TRUE))
  plot(twogenes_df_ordered$Lim1, ylab = "", xlab = "", pch=16, cex = 0.5, main = "Lim1")
  plot(twogenes_df_ordered$gene2, ylab = "", xlab = "", pch=16, cex = 0.5, main = paste(rownames(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[i,])))
  dev.off()
  print(i)
}

#Manual checks:
#1. Lim1 vs. Lim1
#2. Plot another gene manually

###Check Katarina's differentially expd genes
Lawf1_genes <- c("dpr1", "dpr2", "dpr3", "dpr6", "dpr11", "Dscam2", "Fas2", "kek1", "kirre")
Lawf2_genes <- c("beat-IIIb", "parvin", "shakB", "DIP-beta", "dpr8", "dpr13", "dpr17", "side")
Lawf1_gene_rows <- c()
Lawf2_gene_rows <- c()
for (i in 1:length(Lawf1_genes)){
  Lawf1_gene_rows[i] <- which(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly$GeneSymbol_Flybase == paste(Lawf1_genes[i]))}
for (i in 1:length(Lawf2_genes)){
  Lawf2_gene_rows[i] <- which(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly$GeneSymbol_Flybase == paste(Lawf2_genes[i]))}

#plot Katarina's Lawf1 genes
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180601_LawfClustersCloserLook/Other_resulting_Plots/")
quartz("title", 6, 10)
layout(matrix(c(1:11), nrow = 11, byrow = TRUE))
par(mar=c(0,4,1,1), oma = c(0, 0, 2, 0), xpd = NA)
plot(twogenes_df_ordered$Lim1, ylab = "", xlab = "", pch=16, cex = 0.5, xaxt = "n")
legend("topleft", "Lim1",  bty="n", text.col = 'red', text.font = 2, cex = 1.2)
title("Genes >5X enriched in adult Lawf1 neurons", line = 1)
for (i in 1:length(Lawf1_gene_rows))
{
  twogenes_df <- data.frame(as.numeric(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lim1_row, 4:855]), 
                            as.numeric(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lawf1_gene_rows[i], 4:855]))
  colnames(twogenes_df) <- c("Lim1", "gene2")
  #reorder based on Lim1 expression 
  twogenes_df_ordered<-twogenes_df[order(twogenes_df$Lim1),]
  plot(twogenes_df_ordered$gene2, ylab = "", xlab = "", pch=16, cex = 0.5, xaxt = "n")
  legend("topleft", paste(rownames(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lawf1_gene_rows[i],]), "-",
                          k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lawf1_gene_rows[i],3]), bty="n", text.col = 'red', text.font = 2, cex = 1.2) 
}
axis(1)
mtext("Log10(Normalized UMI)", side = 2,  srt=90, line=2, at = c(-90,5), cex = 0.7)
quartz.save(file = "Lawf1Genes.png", type = "png", device = dev.cur()) 
            
#plot Katarina's Lawf2 genes
quartz("title", 6, 10)
layout(matrix(c(1:11), nrow = 11, byrow = TRUE))
par(mar=c(0,4,1,1), oma = c(0, 0, 2, 0), xpd = NA)
plot(twogenes_df_ordered$Lim1, ylab = "", xlab = "", pch=16, cex = 0.5, xaxt = "n")
legend("topleft", "Lim1",  bty="n", text.col = 'red', text.font = 2, cex = 1.2)
title("Genes >5X enriched in adult Lawf2 neurons", line = 1)
for (i in 1:length(Lawf2_gene_rows))
{
  twogenes_df <- data.frame(as.numeric(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lim1_row, 4:855]), 
                            as.numeric(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lawf2_gene_rows[i], 4:855]))
  colnames(twogenes_df) <- c("Lim1", "gene2")
  #reorder based on Lim1 expression 
  twogenes_df_ordered<-twogenes_df[order(twogenes_df$Lim1),]
  plot(twogenes_df_ordered$gene2, ylab = "", xlab = "", pch=16, cex = 0.5, xaxt = "n")
  legend("topleft", paste(rownames(k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lawf2_gene_rows[i],]), "-",
                          k20_Lawf_cluster_GBM_wGeneNames_ExpdGenesOnly[Lawf2_gene_rows[i],3]), bty="n", text.col = 'red', text.font = 2, cex = 1.2) 
}
axis(1)
mtext("Log10(Normalized UMI)", side = 2,  srt=90, line=2, at = c(-90,5), cex = 0.7)
quartz.save(file = "Lawf2Genes.png", type = "png", device = dev.cur()) 

####Pull out expressed TF, cell adhesion molecules, etc.
