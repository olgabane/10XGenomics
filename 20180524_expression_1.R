#20180524_expression_1.R
#Load gene-barcode-matrix two ways - from raw data and using cellranger kit. 
#Verify that matrix elements are equivalent. Convert sparse to full matrix.
#Normalize matrix two ways - using cellranger kit and manually from raw data.
#20180531 update: this code doesn't do exactly what it claims to do. I manually normalize matrix that result from cellranger, but not matrix.mtx file
#For normalization of matrix.mtx file, see "20180531_GBMpulled.R"

#load sparse matrix (.mtx) file
setwd("/Users/Olga/Downloads/cellranger/count-Neuron/outs/filtered_gene_bc_matrices/cellranger/")
library(Matrix)
m1 <- readMM('matrix.mtx')

#Loads matrix.mtx, but with gene names (rows) and cell barcodes (columns) included.
library(cellrangerRkit)
cellranger_pipestance_path <- "/Users/Olga/Downloads/cellranger/count-Neuron"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
m2 <- exprs(gbm)

#convert sparse matrices to full matrices
m3 <- as.matrix(m2)
m4 <- as.matrix(m1)

#verify that matrices are identical, except for column and row names in m3
all(m3 == m4) #comparison comes out true as expected
#sanity check: if i change value of an element in m3, the two matrices are no longer equal. 
m5 <- m3
m5[1,1] <- 1
all(m3 == m5)

#continue with cellRangerR kit to normalize matrix
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))
m1_norm <- exprs(gbm_bcnorm)
m1_log <- exprs(gbm_log)

#note that normaålized matrix keeps only genes expressed in one of 7927 cells 
print(dim(m1))
#result: [1] 17433  7027
print(dim(m1_norm))
#result: [1] 11954  7027
print(dim(m1_log))
#result: [1] 11954  7027

#attempt to understand how to get to normalized matrix
m1_norm_fullmatrix <- as.matrix(m1_norm)
m1_log_fullmatrix <- as.matrix(m1_log)
m1_norm_fullmatrix[8:15, 1:3]
m1_log_fullmatrix[8:15, 1:3]

#based on output of the last two lines, log transformation calculates log10(value+1)
#log transform manually and verify that matrices are equivalent to make sure I understand normalization
all(m1_log_fullmatrix == m1_norm_fullmatrix) #not equal, as expected
m_normlogmanual <-log10(m1_norm_fullmatrix+1)
all(m1_log_fullmatrix == m_normlogmanual) #true as expected, so I understand normalization

#try to normalize original matrix (m3) manually to make sure I understand normalization
#sum UMIs for all cells (store in vector x)
x<-c()
for(i in 1:7027)
  x[i]<-sum(m3[,i])
#get median of UMI values
med_x <- median(x)
m3_norm <- m3
#divide each value in column by sum of UMIs
for(i in 1:7027)
  for(j in 1:17433)
    m3_norm[j,i] <- m3_norm[j,i]/x[i]
#multiply by median UMI count
m3_norm <- med_x*(m3_norm)

#delete all rows(genes) which have values of 0
m3_norm_nozeros <- m3_norm[which(rowSums(m3_norm) > 0),]
all(dim(m3_norm_nozeros) == dim(gbm_log)) #output: TRUE (confirms that deleted correct amt of rows)
#HAVE I SUCCESSFULLY NORMALIZED GENE-BARCODE-MATRIX??
all(m3_norm_nozeros == m1_norm_fullmatrix) #output: FALSE. So matrix still not fully normalized
all.equal(m3_norm_nozeros, m1_norm_fullmatrix) #output: TRUE
which(which(m3_norm_nozeros == m1_norm_fullmatrix) == FALSE) #output: TRUE
#all.equal and which functions claim matrices are identical. 
#but all() does not.
#if matrices identical sum should be 0. BUT it's -6.4506261443497692198e-11
sum(m1_norm_fullmatrix-m3_norm_nozeros)
#to find out why, I found one non-zero value in resulting matrix.
#If I go out 20 digits,there is a difference. So really, the matrices are effectively the same.  
text <- sum(m1_norm_fullmatrix-m3_norm_nozeros)
options(digits=20)
m1_norm_fullmatrix[8,3]
#result: [1] 0.94528228924980661763
m3_norm_nozeros[8, 3]
#result: [1] 0.94528228924980672865

#Log transform normalized matrix (see above for calculation)
GeneBarcodeNormaLogMatrix_Final <-log10(m3_norm_nozeros+1)
all.equal(GeneBarcodeNormaLogMatrix_Final, m1_log_fullmatrix) #get TRUE as expected

#save manually made normalized, log transformed Gene-Barcode-Normalized-Logged Matrix
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180524 Analysis_Gene expression")
write.csv(GeneBarcodeNormaLogMatrix_Final, "GeneBarcodeNormaLogMatrix_Final.csv")

