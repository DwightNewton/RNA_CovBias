#Picard output stats
library(reshape2)
library(ggplot2)
library(multtest)



#plot histograms, and sequencing bias measures (5' vs 3', 5', and 3')
histodata <- read.csv("Picard_human_histogram_output.csv", row.names = "X")


#get group labels
ptable <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/ptable.csv")
ptable_small <- ptable[,c(2,4)]
ptable_small <- ptable_small[!duplicated(ptable_small$`HU.`),]

#Clean up subject ID labels and other inconsistencies 
histodata$ID <- as.numeric(gsub("^.*-Hu+", "", histodata$Library))
histodata <- merge(histodata, ptable_small, by.x="ID", by.y="HU.", all.x=TRUE)

names(histodata)[103] <- "Group"

histodata$CellType <- gsub("-Hu.*$", "", histodata$Library)
unique(histodata$CellType)

histodata$CellType[histodata$CellType=="PV"] <- "PVALB"

histodata$Hybrid <- paste0(histodata$CellType, "_", histodata$Group)


#Melt data for easier handling
histodata_melt <- melt(histodata[,-c(1,2)])
histodata_melt$variable <- as.numeric(gsub("Position_", "", histodata_melt$variable))

#Calculate mean and SD for each subgroup
histo_means <- aggregate(histodata_melt[], by=list(histodata_melt$Hybrid, histodata_melt$variable), FUN="mean")
histo_SDs <- aggregate(histodata_melt[], by=list(histodata_melt$Hybrid, histodata_melt$variable), FUN="sd")

histo_means <- histo_means[,c(1,2,7)]
histo_SDs <- histo_SDs[,c(1,2,7)]

names(histo_means) <- c("Hybrid", "Position", "Coverage")
names(histo_SDs) <- c("Hybrid", "Position", "SD")

histo_summary <- merge(histo_means, histo_SDs)
histo_summary$CellType <- gsub("_.*", "", histo_summary$Hybrid)
histo_summary$Group <- gsub(".*_", "", histo_summary$Hybrid)

histo_summary$Group <- factor(histo_summary$Group, levels = c("Control", "MDD", "Bipolar", "SCHIZ"))
histo_summary$CellType <- factor(histo_summary$CellType, levels = c("Pyr-L2n3", "Pyr-L5n6", "SST", "PVALB", "VIP"))


histo_summary_CTAVG <-aggregate(histo_summary, by=list(histo_summary$CellType, histo_summary$Position), FUN="mean")
histo_summary_CTAVG$CellType <- gsub("_.*", "", histo_summary_CTAVG$Group.1)
histo_summary_CTAVG$Group <- gsub(".*_", "", histo_summary_CTAVG$Group.1)


#Plotting individual rather than facets for presentation purposes (too crowded for facets)
#Histograms by cell-type: controls
Con <- ggplot(histo_summary[histo_summary$Group == "Control",], aes(x=Position, fill=CellType)) +
  geom_ribbon(aes(ymin=Coverage-SD, ymax=Coverage+SD), alpha=0.3) +
  scale_fill_manual(values=c("#F8766D", "#b01005","#00BFC4","#7CAE00","#C77CFF")) +
  scale_colour_manual(values=c("#F8766D", "#b01005","#00BFC4","#7CAE00","#C77CFF")) +
  geom_line(aes(y=Coverage, colour=CellType), size=1) +
  theme_classic() +
  #scale_y_continuous(limits = c(0.9,1.3)) +
  ggtitle("Coverage Histogram - Controls") +
  theme(legend.position = "none")



#MDD
MDD <- ggplot(histo_summary[histo_summary$Group == "MDD",], aes(x=Position, fill=CellType)) +
  geom_ribbon(aes(ymin=Coverage-SD, ymax=Coverage+SD), alpha=0.3) +
  scale_fill_manual(values=c("#F8766D", "#b01005","#00BFC4","#7CAE00","#C77CFF")) +
  scale_colour_manual(values=c("#F8766D", "#b01005","#00BFC4","#7CAE00","#C77CFF")) +
  geom_line(aes(y=Coverage, colour=CellType), size=1) +
  theme_classic() +
  #scale_y_continuous(limits = c(0.9,1.3)) +
  ggtitle("Coverage Histogram - MDD") +
  theme(legend.position = "none")

#Bipolar
Bipolar <- ggplot(histo_summary[histo_summary$Group == "Bipolar",], aes(x=Position, fill=CellType)) +
  geom_ribbon(aes(ymin=Coverage-SD, ymax=Coverage+SD), alpha=0.3) +
  scale_fill_manual(values=c("#F8766D", "#b01005","#00BFC4","#7CAE00","#C77CFF")) +
  scale_colour_manual(values=c("#F8766D", "#b01005","#00BFC4","#7CAE00","#C77CFF")) +
  geom_line(aes(y=Coverage, colour=CellType), size=1) +
  theme_classic() +
  #scale_y_continuous(limits = c(0.9,1.3)) +
  ggtitle("Coverage Histogram - BPD") +
  theme(legend.position = "none")

#SCHIZ
SCHIZ <- ggplot(histo_summary[histo_summary$Group == "SCHIZ",], aes(x=Position, fill=CellType)) +
  geom_ribbon(aes(ymin=Coverage-SD, ymax=Coverage+SD), alpha=0.3) +
  scale_fill_manual(values=c("#F8766D", "#b01005","#00BFC4","#7CAE00","#C77CFF")) +
  scale_colour_manual(values=c("#F8766D", "#b01005","#00BFC4","#7CAE00","#C77CFF")) +
  geom_line(aes(y=Coverage, colour=CellType), size=1) +
  theme_classic() +
  #scale_y_continuous(limits = c(0.9,1.3)) +
  ggtitle("Coverage Histogram - SCZ") +
  theme(legend.position = "none")

library(gridExtra)
grid.arrange(Con, MDD, Bipolar, SCHIZ, nrow=1)


#Compare distributions across cell-types, disorders, any other A/B grouping using KS tests 
#Follow-up significant tests with bonferroni-corrected t-tests of 5' and 3' regions

#KS, t-test, and Rank-Sum tests of distribution similarity
histo_stats <- read.csv("Picard_human_histogram_output.csv", row.names="X")
ptable <- read.csv("C:/Users/dwigh/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/NovaSeq/ptable.csv")
ptable_small <- ptable[,c(2,4)]
ptable_small <- ptable_small[!duplicated(ptable_small$`HU.`),]

histo_stats$ID <- as.numeric(gsub("^.*-Hu+", "", histo_stats$Library))
histo_stats <- merge(histo_stats, ptable_small, by.x="ID", by.y="HU.", all.x=TRUE)

names(histo_stats)[103] <- "Group"

histo_stats$CellType <- gsub("-Hu.*$", "", histo_stats$Library)
unique(histo_stats$CellType)

histo_stats$CellType[histo_stats$CellType=="PV"] <- "PVALB"

histo_stats$Hybrid <- paste0(histo_stats$CellType, "_", histo_stats$Group)


histo_KS <- aggregate.data.frame(histo_stats, by=list(histo_stats$Hybrid), FUN="mean")
histo_KS <- histo_KS[,-c(2,3)]
histo_KS$CellType <- gsub("_.*", "", histo_KS$Group.1)
histo_KS$Group <- gsub(".*_", "", histo_KS$Group.1)

pairlist <- combn(histo_KS$Group.1, 2)

#Only contrasts of interest
pairlist <- pairlist[,(gsub(".*_","",pairlist[1,]) == gsub(".*_","",pairlist[2,])) |( gsub("_.*","",pairlist[1,]) == gsub("_.*","",pairlist[2,]))]

KS_results <- as.data.frame(matrix(ncol=3, nrow=ncol(pairlist)))
names(KS_results) <- c("Group_1", "Group_2", "p-value")

for(i in 1:ncol(pairlist)){
  group1 <- pairlist[1,i]
  group2 <- pairlist[2,i]
  print(paste(group1, "v.s.", group2, "KS test"))
  KS_res <- ks.test(as.numeric(histo_KS[histo_KS$Group.1 == group1, 4:103]), as.numeric(histo_KS[histo_KS$Group.1 == group2, 4:103]))
  KS_results$Group_1[i] <- group1
  KS_results$Group_2[i] <- group2
  KS_results$`p-value`[i] <- KS_res$p.value
  print(KS_res)
  
}

#T-test at each position for 5'-most 25% and 3'-most 25%
histo_ttests <- histo_stats


ttest_results <- as.data.frame(matrix(ncol=102, nrow=70))
names(ttest_results) <- c("Group_1", "Group_2", paste0("Position_", c(1:100)))


for(i in 1:ncol(pairlist)){
  group1 <- pairlist[1,i]
  group2 <- pairlist[2,i]
  print(paste(group1, "v.s.", group2, "t-test"))
  ttest_results[i, c(1:2)] <- c(group1, group2)
  for(j in 1:100){
    t_res <- t.test(histo_ttests[histo_ttests$Hybrid == group1, j+2], histo_ttests[histo_ttests$Hybrid == group2, j+2])
    ttest_results[i,j+2] <- t_res$p.value
  }
}


#Bonferroni for stringency
ttest_results_bonferroni <- ttest_results

for(i in 1:nrow(ttest_results_bonferroni)){
  ttest_results_bonferroni[i,3:102] <- p.adjust(ttest_results_bonferroni[i,3:102], method = "BH")
}

#Mean difference for each comparison
ttest_means <- histo_stats
ttest_means <- aggregate(ttest_means, by=list(histo_stats$Hybrid), FUN="mean")
ttest_means$CellType <- gsub("_.*", "", ttest_means$Group.1)
ttest_means$Group <- gsub(".*_", "", ttest_means$Group.1)
ttest_means <- ttest_means[,-2]

mean_outputs <- as.data.frame(matrix(ncol=102, nrow=70))
names(mean_outputs) <- c("Group_1", "Group_2", paste0("Position_", c(11:30, 71:90)))
for(i in 1:ncol(pairlist)){
  group1 <- pairlist[1,i]
  group2 <- pairlist[2,i]
  print(paste(group1, "v.s.", group2, "t-test"))
  mean_outputs[i, c(1:2)] <- c(group1, group2)
  for(j in 1:100){
    mean_res <- ttest_means[ttest_means$Group.1 == group1, j+2] - ttest_means[ttest_means$Group.1 == group2, j+2]
    mean_outputs[i,j+2] <- mean_res
  }
}

write.csv(mean_outputs, "Coverage ttest MEANS.csv")
write.csv(ttest_results_bonferroni, "Bonferroni-corrected coverage ttests.csv")
write.csv(ttest_results, "Coverage ttests.csv")


