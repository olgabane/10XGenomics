#20180721_oldneuronscompare_GOterms.R

setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180706_oldneuronscompare")

#Got GO terms for relative_expresssion_sorted_OMannotated FBgn IDs from flybase: http://flybase.org/batchdownload. Choose "Precomputed files" from data source!
#Used FBgn IDs for 5X-enriched genes;  
#There are 204 Lim1- cells, 292 Lim1+ cells. Used only genes expressed in > 10% of cells.
a<-read.csv("Enriched_Lim1neg_Goterms.txt", sep = "\t", header = FALSE)
Enriched_Lim1neg_Goterms <- data.frame(a$V2, a$V5)
colnames(Enriched_Lim1neg_Goterms) <-c("FlyBaseID", "GO_Term")

a<-read.csv("Enriched_Lim1pos_Goterms.txt", sep = "\t", header = FALSE)
Enriched_Lim1pos_Goterms <- data.frame(a$V2, a$V5)
colnames(Enriched_Lim1pos_Goterms) <-c("FlyBaseID", "GO_Term")

# load the GO library, get list of GO ID decriptions
source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
library(GO.db)
goterms <- Term(GOTERM)
goterms <- as.data.frame(goterms)
goterms$GO_ID <- rownames(goterms)
rownames(goterms) <- c()

#Add goterm description to FBgn IDs enriched in Lim1pos, Lim1neg
pos_rows <- match(Enriched_Lim1pos_Goterms$GO_Term, goterms$GO_ID)
neg_rows <- match(Enriched_Lim1neg_Goterms$GO_Term, goterms$GO_ID)

Enriched_Lim1pos_Goterms_defined <- data.frame(Enriched_Lim1pos_Goterms, goterms[pos_rows,])
Enriched_Lim1neg_Goterms_defined <- data.frame(Enriched_Lim1neg_Goterms, goterms[neg_rows,])

#Save files
write.csv(Enriched_Lim1pos_Goterms_defined, "Enriched_Lim1pos_Goterms_defined.csv")
write.csv(Enriched_Lim1neg_Goterms_defined, "Enriched_Lim1neg_Goterms_defined.csv")


