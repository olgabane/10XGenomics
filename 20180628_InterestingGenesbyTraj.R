#20180628_InterestingGenesbyTraj.R
#June 28, 2018
#Manually encoded trajectory over time for all 8032 genes
#This script: pull out genes with certain trajectory patterns

#load manually encoded trajectory as data frame
setwd("/Users/Olga/Google Drive/Desplan Lab/Notebooks/Notebook 5/10X processing/20180619_Monocle/PC3/AllGenesOverPseudotime/")
df_AllGenesOverPseudotime<-read.csv("AllGenesOverPseudotime_Summary_R.csv")
#clean up df
df_AllGenesOverPseudotime$X <- df_AllGenesOverPseudotime$X.1 <- df_AllGenesOverPseudotime$X.2 <-NULL
df_AllGenesOverPseudotime[is.na(df_AllGenesOverPseudotime)] <- 0

#subset df to include only genes that change over time in neurons
df_ChangedGenesOverPseudotime <- df_AllGenesOverPseudotime[which(df_AllGenesOverPseudotime$Neuron_Down == 1 | 
                                                                   df_AllGenesOverPseudotime$Neuron_Up == 1),]
#Bin trajectories of changed genes into descriptive categories
Trajectory_category<-function(a, b, c){
  which(df_ChangedGenesOverPseudotime$Neuron_Straight == a & 
                     df_ChangedGenesOverPseudotime$Neuron_Up == b &
                     df_ChangedGenesOverPseudotime$Neuron_Down == c)  
}
Weak_up <- Trajectory_category(1, 1, 0)
Weak_down <- Trajectory_category(1, 0, 1)
Weak_up_down <- Trajectory_category(1, 1, 1)
Up <- Trajectory_category(0, 1, 0)
Down <- Trajectory_category(0, 0, 1)
Up_Down <- Trajectory_category(0, 1, 1)
#verify that binning gives expected number of instances
dim(df_ChangedGenesOverPseudotime)[1] == sum(length(Weak_up), length(Weak_down), length(Weak_up_down), length(Up), length(Down), length(Up_Down))
#Add descriptive catergories to df
df_ChangedGenesOverPseudotime[Weak_up, 11] <- "Weak_up"
df_ChangedGenesOverPseudotime[Weak_down, 11] <- "Weak_down"
df_ChangedGenesOverPseudotime[Weak_up_down, 11] <- "Weak_up_down"
df_ChangedGenesOverPseudotime[Up, 11] <- "Up"
df_ChangedGenesOverPseudotime[Down, 11] <- "Down"
df_ChangedGenesOverPseudotime[Up_Down, 11] <- "Up_Down"
colnames(df_ChangedGenesOverPseudotime)[11] <- "Trajectory"

#Identify TFs, chromatin related protein, and transcriptional machinery components among changed genes
#Rheet et al., 2015, Cell Reports, TableS1: 
#"We surveyed the literature and gathered a list of 996 genes, containing TFs with characterized binding domains, 
#computationally predicted (putative) TFs, chromatin-related proteins and transcriptional machinery components"
TFList <- read.csv("Rhee2015_TableS1.csv")
TFList[, 3:7] <- NULL
match<-match(df_ChangedGenesOverPseudotime$Submitted_ID, TFList$Flybase.ID)
match[is.na(match)] <- 0
match[which(match > 0)] <- 1
df_ChangedGenesOverPseudotime$TF_Rhee2015 <- match

#Add automated gene summaries from flybase to resulting changed genes
GeneSummaries_FB <- read.csv("automated_gene_summaries.tsv", sep = "\t", header = FALSE)
colnames(GeneSummaries_FB) <- c("Flybase_ID", "Gene_summary")
match_FBID<-match(GeneSummaries_FB$Flybase_ID, df_ChangedGenesOverPseudotime$Submitted_ID, nomatch = 0)
match_FBID[which(match_FBID > 0)] <- 1
GeneSummaries_FB_subset <- GeneSummaries_FB[which(match_FBID == 1),]
#re-order, clean up to be able to add to df_ChangedGenesOverPseudotime
for (i in 1:13){
GeneSummaries_FB_subset[nrow(GeneSummaries_FB_subset)+1,] <- NA
i= i+1
}
GeneSummaries_FB_subset$IDs_in_f_ChangedGenesOverPseudotime <- df_ChangedGenesOverPseudotime$Submitted_ID
rownames(GeneSummaries_FB_subset) <- 1:nrow(GeneSummaries_FB_subset)
a<-match(GeneSummaries_FB_subset$IDs_in_f_ChangedGenesOverPseudotime, GeneSummaries_FB_subset$Flybase_ID)
nomatch<-which(is.na(a))

#clean up: insert column for each ID that has no corresponding summary in GeneSummaries_FB
df1 <- GeneSummaries_FB_subset[1:(nomatch[1]-1),]
df1[nrow(df1)+1, ] <- NA

for(i in 2:length(nomatch)){
assign(paste("df", i, sep = ""), GeneSummaries_FB_subset[(nomatch[i-1] - (i-2)):(nomatch[i]-1-(i-1)),])
}

dfList <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13)

for(i in 2:length(nomatch)){
dfList[[i]][nrow(dfList[[i]])+1, ] <- NA
assign(paste("df", i, sep = ""), dfList[[i]])
}

df14 <- GeneSummaries_FB_subset[(nomatch[13]-(12)):(dim(GeneSummaries_FB_subset)[1]),]

GeneSummaries_FB_clean <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14)

rownames(GeneSummaries_FB_clean) <- 1:nrow(GeneSummaries_FB_clean)
GeneSummaries_FB_clean<- GeneSummaries_FB_clean[-(1577:1589), ]

#re-append ID's from df_ChangedGenesOverPseudotime to verify that summaries are matched correctly
GeneSummaries_FB_clean$IDs_in_f_ChangedGenesOverPseudotime <- df_ChangedGenesOverPseudotime$Submitted_ID

#once verified, append gene summaries to df_ChangedGenesOverPseudotime
df_ChangedGenesOverPseudotime$Gene_summary <- GeneSummaries_FB_clean$Gene_summary

#Save resulting df
write.csv(df_ChangedGenesOverPseudotime, "AllGenesOverPseudotime_ProcessedJune29.csv")
