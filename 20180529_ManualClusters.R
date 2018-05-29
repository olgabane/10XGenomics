#tsne

library(gatepoints)

setwd("/Users/Olga/Downloads/cellranger/count-Neuron/outs/analysis/tsne/2_components/")
m<-read.csv("projection.csv")

#change dir to where output from fhs() will be saved
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180529_ClusterIDGeneExpr/Dataframes_from_gatepoints/")

#Draw plot, choose points
#Convert rows containing chosen points to numerics, select only those points for tsne plot 
#save chosen points because selection is "fleeting" (although bounds of selection are saved)
#replot with chosen points
quartz("tsne", 5, 5.3)
plot(m$TSNE.1, m$TSNE.2, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.3, bg = 'black')
as.numeric(cells_selected)<-fhs(m[, 2:3], mark = FALSE, names = TRUE)
write.csv(cells_selected, "cells_selected_repo.csv")
m_sel <- m[cells_selected,]
quartz("tsne", 5, 5.3)
plot(m$TSNE.1, m$TSNE.2, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.3, bg = 'black')
points(m_sel$TSNE.1, m_sel$TSNE.2,  pch = 16, cex = 0.3, col = 'red')

#Run above script, re-drawing 