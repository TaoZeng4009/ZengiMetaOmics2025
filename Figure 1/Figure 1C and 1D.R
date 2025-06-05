### Pipeline for standard analyses of alpha diversity of bacterial communities
### Based on R v4.0.2
library(vegan)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ape)
library(picante)
library(dplyr)
library(multcomp)

### Loading Parameters & Processing
source("Functions/Dataprocess/Loading.R")

cbbPalette <- Loading_Color()
bac <- Loading_Table("Input/taxa_table.xls")
group <- Loading_Group("Input/group.txt")
tree <- Loading_Tree("Input/phylogeny.tre")

result <- Data_Process(bac,group,tree)
bac1 <- result[[1]]
bac2 <- result[[2]]
bac3 <- result[[3]]
tree <- result[[4]]
Group_numb <- result[[5]]
Sample_numb <- result[[6]]
wid <- result[[7]]

## Alpha diversity analysis
dir.create("Results")
dir.create("Results/02.Alpha")

### Alpha diversity indices calculation
dir.create("Results/02.Alpha/01.Estimators")
source("Functions/Alpha/alpha_diversity_indices.R")
alpha.bac <- Alpha_diversity_index(bac1,tree)
alpha <- alpha.bac[,-c(3,5,9)]
colnames(alpha) <- c("Observed_species","Chao1","ACE","Shannon","Simpson","Pielou_J",
                     "Pd_faith")
alpha$ACE[is.na(alpha$ACE)] <- alpha$Chao1[is.na(alpha$ACE)]
alpha1 <- data.frame(ID = rownames(alpha),alpha)

write.table(alpha1,"Results/02.Alpha/01.Estimators/alpha_diversity_indices.txt",
            sep = "\t",row.names = FALSE)

### Differences of alpha diversity indices among different groups
#### Boxplot
dir.create("Results/02.Alpha/01.Estimators/alpha_div_boxplot")
#### Jitter
dir.create("Results/02.Alpha/01.Estimators/alpha_div_jitter")
source("Functions/Alpha/Diff_alpha2.R")
result <- diff.alpha.jitter(alpha1,group)
for (i in 2:ncol(alpha1)) {
    pdf(file = paste("Results/02.Alpha/01.Estimators/alpha_div_jitter/",
                     colnames(alpha1)[i],".pdf",sep = ""),
        width = 0.8 + 0.6*Group_numb,height = 3.6)
    print(result[[i-1]])
    dev.off()
}


