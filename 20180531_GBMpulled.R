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
#note: CANNOT use all() or all.equal() due to different attributes. Must compare values
a<-which(m1 != m2)
b<-which(m1 == m2)
c<-which(round(m1, digits = 10) != round(m2, digits = 10))
d<-which(round(m1, digits = 10) == round(m2, digits = 10))
sprintf("a = %i, b = %i, c = %i, d = %i, a+b = %i, c+d= %i", length(a), length(b), length(c), length(d), length(a)+length(b), length(c)+length(d))
sprintf("c+d = %i %i", length(c)+length(d))
#second way to verify that matrices are equal.
sum(round(m1, digits = 10) - round(m2, digits = 10))



