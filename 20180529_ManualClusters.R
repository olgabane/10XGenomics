#tsne

library(gatepoints)

setwd("/Users/Olga/Downloads/cellranger/count-Neuron/outs/analysis/tsne/2_components/")
m<-read.csv("projection.csv")

#change dir to where output from fhs() will be saved
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180529_ClusterIDGeneExpr/Dataframes_from_gatepoints/")

#Function: Draw plot, choose points, convert rows corresp to chosen pts to numeric, save chosen pts (rows)
open_tsne_plot_for_pt_selection <- function(x)
{
  quartz("tsne_original", 5, 5.3)
  plot(m$TSNE.1, m$TSNE.2, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.3, bg = 'black')
  cells_selected <- fhs(m[, 2:3], mark = FALSE, names = TRUE)
  cells_selected <- as.numeric(cells_selected)
  write.csv(cells_selected, file = x)
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

#Draw repo cluster
open_tsne_plot_for_pt_selection("Repo_rows.csv")
plot_selection_over_tsne("Repo_cluster_plot.jpg")







