#20180601_LawfClusterCloserLook.R
#June 1, 2018, Continued June 4, 2018

#Load GBM file created 20180531 (see 20180531_GBMpulled.R)
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180531_GBMpulled/")
k20_Lawf_cluster_GBM_wGeneNames<-read.table("k20_Lawf_cluster_GBM_wGeneNames.txt")

#Plot Normalized UMI values 
gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "Lim1")
quartz("title", 10, 10)
layout(matrix(c(1:8), nrow = 4, byrow = TRUE), widths=c(4,1))
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), col=ifelse((as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]) == 0), "red", "blue"), pch = 16, cex = 0.5, main="Lim1", xlab = "852 'Lawf' cells", ylab = "Log10(Normalized UMI)")
p2=barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0))), main = "Lim1", col = c("red", "blue"))
axis(1, at = p2, labels = c("0", ">0"))

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "hth")
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), col=ifelse((as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]) == 0), "red", "blue"), pch = 16, cex = 0.5, main="hth", xlab = "852 'Lawf' cells", ylab = "Log10(Normalized UMI)")
p2=barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0))), main = "hth", col = c("red", "blue"))
axis(1, at = p2, labels = c("0", ">0"))

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "hbn")
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), col=ifelse((as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]) == 0), "red", "blue"), pch = 16, cex = 0.5, main="hbn", xlab = "852 'Lawf' cells", ylab = "Log10(Normalized UMI)")
p2=barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0))), main = "hbn", col = c("red", "blue"))
axis(1, at = p2, labels = c("0", ">0"))

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "eya")
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), col=ifelse((as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]) == 0), "red", "blue"), pch = 16, cex = 0.5, main="eya", xlab = "852 'Lawf' cells", ylab = "Log10(Normalized UMI)")
p2=barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0))), main = "eya", col = c("red", "blue"))
axis(1, at = p2, labels = c("0", ">0"))

###Plot above data as heat map to look at expression of all 4 factors in single cell
Lim1_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "Lim1")
hth_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "hth")
hbn_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "hbn")
eya_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "eya")
fourgenes_df <- data.frame(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[Lim1_row, 4:855]), as.numeric(k20_Lawf_cluster_GBM_wGeneNames[hth_row, 4:855]), as.numeric(k20_Lawf_cluster_GBM_wGeneNames[hbn_row, 4:855]), as.numeric(k20_Lawf_cluster_GBM_wGeneNames[eya_row, 4:855]))
colnames(fourgenes_df) <- c("Lim1", "hth", "hbn", "eya")
#reorder based on Lim1 expression 
fourgenes_df_ordered<-fourgenes_df[order(fourgenes_df$Lim1),]
#set bins for cut function
bins <- seq(range(fourgenes_df_ordered)[1], range(fourgenes_df_ordered)[2], by = range(fourgenes_df_ordered)[2]/7)
bins[1] <- -1 #because 0 won't be inclusive
Lim1_colors<-rainbow(7, start = 0, end = 2/3)[cut(fourgenes_df_ordered$Lim1, breaks = bins)]
hth_colors<-rainbow(7, start = 0, end = 2/3)[cut(fourgenes_df_ordered$hth, breaks = bins)]
hbn_colors<-rainbow(7, start = 0, end = 2/3)[cut(fourgenes_df_ordered$hbn, breaks = bins)]
eya_colors<-rainbow(7, start = 0, end = 2/3)[cut(fourgenes_df_ordered$eya, breaks = bins)]

quartz("title", 10, 3)
layout(matrix(c(1:8), nrow = 4, byrow = TRUE), width = c(4,2))
par(mar = c(0,2,0,0), oma =  c(6,0,1,1), xpd = NA)
barplot(rep(1, 852), col = Lim1_colors, border =Lim1_colors, xlab = "852 'Lawf' cells", axes = FALSE)
title(ylab = "Lim1", line = -1, cex.lab = 1.5)
plot(fourgenes_df_ordered$Lim1, col=Lim1_colors, ylab = "", xlab = "", pch=16, xaxt = 'n', cex = 0.5, axes = FALSE)
box(lty=1)
axis(2, labels = c("0", "1"), at = c(0, 1))
barplot(rep(1, 852), col = hth_colors, border =hth_colors, xlab = "852 'Lawf' cells", axes = FALSE)
title(ylab = "hth", line = -1, cex.lab = 1.5)
plot(fourgenes_df_ordered$hth, col=hth_colors, ylab = "", xlab = "", pch=16, xaxt = 'n', cex = 0.5, axes = FALSE)
box(lty=1)
axis(2, labels = c("0", "1"), at = c(0, 1))
barplot(rep(1, 852), col = hbn_colors, border =hbn_colors, xlab = "852 'Lawf' cells", axes = FALSE)
title(ylab = "hbn", line = -1, cex.lab = 1.5)
plot(fourgenes_df_ordered$hbn, col=hbn_colors, ylab = "", xlab = "", pch=16, xaxt = 'n', cex = 0.5, axes = FALSE)
box(lty=1)
axis(2, labels = c("0", "1"), at = c(0, 1))
title(ylab = "Log10(Normalized UMI)", line = 2, adj = 0.2)
barplot(rep(1, 852), col = eya_colors, border =eya_colors, axes = FALSE)
title(ylab = "eya", line = -1, cex.lab = 1.5)
title(xlab = "852 'Lawf' cells")
plot(fourgenes_df_ordered$eya, col=eya_colors, xlab = "852 'Lawf' cells", ylab = "", pch=16, cex = 0.5, axes = FALSE)
box(lty=1)
axis(1, labels = c("0","852"), at = c(0, 852))
axis(2, labels = c("0", "1"), at = c(0, 1))

###Plot above data a third way: most compact, gives most info
quartz("title", 6, 4)
par(mar = c(0,4,0,0), oma =  c(6,0,2,1), xpd = NA)
layout(matrix(c(1, 5, 
                2, 6, 
                3, 7, 
                4, 8), nrow = 4, byrow = TRUE), widths = c(2.5, 1))
plot(fourgenes_df_ordered$Lim1, col=ifelse(fourgenes_df_ordered$Lim1 == 0, "red", "blue"), ylab = "", xlab = "", pch=16, xaxt = 'n', cex = 0.5, ylim=c(0,1.8), yaxt='n')
axis(2, at = c(0, 1.8), labels = c("0", "1.8"), las = 1)
plot(fourgenes_df_ordered$hth, col=ifelse(fourgenes_df_ordered$hth == 0, "red", "blue"), ylab = "", xlab = "", pch=16, xaxt = 'n', cex = 0.5, ylim=c(0,1.8), yaxt='n')
plot(fourgenes_df_ordered$hbn, col=ifelse(fourgenes_df_ordered$hbn == 0, "red", "blue"), ylab = "", xlab = "", pch=16, xaxt = 'n', cex = 0.5, ylim=c(0,1.8), yaxt='n')
plot(fourgenes_df_ordered$eya, col=ifelse(fourgenes_df_ordered$eya == 0, "red", "blue"), ylab = "", xlab = "", pch=16, xaxt = 'n', cex = 0.5, ylim=c(0,1.8), yaxt='n')
axis(1, at = c(1,  852), labels=c("0", "852"))
title(xlab = "852 'Lawf' cells", line = 2)

axisbox <- function()
{
  axis(2, at = c(0,  800), labels=c("0", "800"), las = 1)
  box(lty = 1)
}

Lim1_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "Lim1")
hth_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "hth")
hbn_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "hbn")
eya_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "eya")
#par(mar = c(0,2,0,0), oma =  c(2,0,1,1), xpd = NA,  new = TRUE)
barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[Lim1_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[Lim1_row,] != 0))), 
        col = c("red", "blue"), ylim = c(0, 900), axes = FALSE)
axisbox()
barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[hth_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[hth_row,] != 0))), 
        col = c("red", "blue"), ylim = c(0, 900), axes = FALSE)
box(lty = 1)
barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[hbn_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[hbn_row,] != 0))), 
        col = c("red", "blue"), ylim = c(0, 900), axes = FALSE)
box(lty = 1)
barplot(c(length(which(k20_Lawf_cluster_GBM_wGeneNames[eya_row,] == 0)), length(which(k20_Lawf_cluster_GBM_wGeneNames[eya_row,] != 0))), 
        col = c("red", "blue"), ylim = c(0, 900), axes = FALSE)
title(xlab = "Expression", line = 2)
box(lty = 1)
axis(1, at= c(0.7, 1.9), labels = c("0", ">0"), lty=0)









####NEXT: Delete genes with 0 when summed across all cells
  #how many are left?
#####Sum UMIs of all genes. 
  #Look at most highly exp'd
#####Look at most variable genes across Lawfs. 
  #Is Lim1 one of them?
  #If yes, do any correlate w Lim1?



##############

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "elav")
quartz("title", 10, 10)
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main)

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "repo")
quartz("title", 10, 10)
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main)
a<-which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)
b<-which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0)
length(a)
length(b)
#### delete these if they don't have elav in them!!

gene_row <- which(k20_Lawf_cluster_GBM_wGeneNames$GeneSymbol_Flybase == "ase")
quartz("title", 10, 10)
plot(as.numeric(k20_Lawf_cluster_GBM_wGeneNames[gene_row, 4:855]), pch = 16, cex = 0.5, main = "Lim1 in Lawf cluster")
a<-which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] == 0)
b<-which(k20_Lawf_cluster_GBM_wGeneNames[gene_row,] != 0)