#Manually plot clusters based on known cell markers: ID glia, progenitors and Lawfs

library(gatepoints)

#Function: Draw plot, choose points, convert rows corresp to chosen pts to numeric, save chosen pts (rows)
open_tsne_plot_for_pt_selection <- function(x)
{
  quartz("tsne_original", 5, 5.3)
  plot(m$TSNE.1, m$TSNE.2, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.3, bg = 'black')
  cells_selected <- fhs(m[, 2:3], mark = FALSE, names = TRUE)
  cells_selected <- as.numeric(cells_selected)
  write.csv(cells_selected, file = x)
  return(cells_selected)
}

#Function: plot chosen points over tsne plot
plot_selection_over_tsne <- function(x)
{
  m_sel <- m[cells_selected,]
  quartz("tsne_highlight", 5, 5.3)
  plot(m$TSNE.1, m$TSNE.2, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.3, bg = 'black')
  points(m_sel$TSNE.1, m_sel$TSNE.2,  pch = 16, cex = 0.3, col = 'red')
  dev.copy(png, file = x)
  dev.off()
}


setwd("/Users/Olga/Downloads/cellranger/count-Neuron/outs/analysis/tsne/2_components/")
m<-read.csv("projection.csv")

#change dir to where output from fhs() will be saved
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180529_ClusterIDGeneExpr/Dataframes_from_gatepoints/")

#Draw repo cluster
cells_selected <- open_tsne_plot_for_pt_selection("Repo_rows.csv")
plot_selection_over_tsne("Repo_cluster_plot.png")

#Draw progenitor cluster
cells_selected <- open_tsne_plot_for_pt_selection("Progenitor_rows.csv")
plot_selection_over_tsne("Progenitor_cluster_plot.png")

#Draw Lawf cluster
cells_selected <- open_tsne_plot_for_pt_selection("Lawf_rows.csv")
plot_selection_over_tsne("Lawf_cluster_plot.png")

#Draw overlap of repo and progenitor cluster
repo_df <- read.csv("Repo_rows.csv")
progenitor_df <- read.csv("Progenitor_rows.csv")
intersection <- intersect(repo_df$x, progenitor_df$x)
lawf_df <- read.csv("Lawf_rows.csv")
m_sel1 <- m[repo_df$x,]
m_sel2 <- m[progenitor_df$x,]
m_sel3 <- m[intersection,]
m_sel4 <- m[lawf_df$x, ]
quartz("tsne_highlight", 5, 5.3)
plot(m$TSNE.1, m$TSNE.2, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.3, bg = 'black')
points(m_sel1$TSNE.1, m_sel1$TSNE.2,  pch = 16, cex = 0.3, col = 'red')
points(m_sel2$TSNE.1, m_sel2$TSNE.2,  pch = 16, cex = 0.3, col = 'blue')
points(m_sel3$TSNE.1, m_sel3$TSNE.2,  pch = 16, cex = 0.3, col = 'purple')
points(m_sel4$TSNE.1, m_sel4$TSNE.2,  pch = 16, cex = 0.3, col = 'green')
text(-49,49, "repo+ - glia", col='red', adj = c(0,0))
text(-49,42, "repo+ progenitor-like", col='purple', adj = c(0,0))
text(-49,35, "progenitor", col='blue', adj = c(0,0))
text(-49,28, "Lawfs", col='green', adj = c(0,0))
dev.copy(png, "repo-prog-Lawf_intersect.png")
dev.off()



