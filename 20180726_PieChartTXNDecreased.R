#Pie chart downregulated genes
#20180726_PieChart_TXNDecreased.R
#For poster
#July 26, 2018

a<-c(8,5,3,20,13,11,23,9,41,8,72,43,62,4,16,104,45,58)
b<-c("ATP", "Calcium binding", "Cell adhesion", "Chromatin regulation", "Cytoskeleton and ECM", "GTPase and GPCR"
,"Kinases and phosphatases"
,"lncRNA"
,"Metabolic enzyme"
,"Mitochondrial proteins"
,"Ribosomal protein"
,"RNA binding and regulation"
,"Transcriptional regulator"
,"Synaptic protein"
,"Transmembrane and cell surface"
,"Unknown"
,"Protein regulator"
,"Other")

df<- data.frame(b, a)

df <- df[order(a),]

quartz("a", 6, 6)
pie(df$a, df$b, col=rainbow(length(df$b)), cex = 0.7, radius = 0.7)
