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
match<-match(df_ChangedGenesOverPseudotime$Submitted_ID, TFList$Flybase.ID )
match[is.na(match)] <- 0
match[which(match > 0)] <- 1
df_ChangedGenesOverPseudotime$TF_Rhee2015 <- match
