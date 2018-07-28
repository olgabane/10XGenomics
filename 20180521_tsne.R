setwd("/Users/Olga/Downloads/cellranger/count-Neuron/outs/analysis/tsne/2_components/") 
m <- read.csv("projection.csv")
x = m$TSNE.1
y = m$TSNE.2
quartz("tsne", 5, 5)
setwd("/Users/Olga/Desktop/Demo/20180521/")
png("tsne.png")
plot(x, y, xlab = "TSNE.1", ylab = "TSNE.2", pch = 16, cex = 0.3, bg = 'black')
dev.off()
