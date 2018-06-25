#20180625_sortplottedgenes.R
#June 25, 2018
#Get gene symbols, names for 8032 genes plotted as expression over pseudotime in 20180621_Monocle_4pseudotrajectories.R script

#Put gene names and symbols corresponding to FBID of plotted genes in dataframe. 
setwd("/Users/Olga/Downloads/")
Converted_FBID <- read.csv("FlyBase_Fields_download.csv")
colnames(Converted_FBID) <- c("Submitted_ID", "Gene_Name", "Gene_Symbol")
FBID_plotted<-read.table("FBID_plotted.txt")
colnames(FBID_plotted) <- c("FBID_plotted")
df<-data.frame(FBID_plotted$FBID_plotted, Converted_FBID$Submitted_ID, Converted_FBID$Gene_Name, Converted_FBID$Gene_Symbol)
colnames(df)<-c("FBID_plotted", "Submitted_ID", "Gene_Name", "Gene_Symbol")
#Sort based on order of saved plots. Get order, sort by match vector, remove column 1
match_vector<-match(df[,2], df[,1])
df$matchvector <- match_vector
df_sorted <- df[order(match_vector),] 
df_sorted[1] <-NULL
write.csv(df_sorted, "df_sorted.csv")
