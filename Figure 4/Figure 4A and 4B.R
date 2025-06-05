### Pipeline for standard analyses of alpha diversity of bacterial communities
### Based on R v4.0.2
library(picante)
library(dplyr)
library(tidyverse)
library(reshape2)
library(linkET)
library(ggraph)
library(tidygraph)
library(igraph)
library(ggplot2)
library(ggrepel)
library(plotly)
library(multcomp)

### Loading Parameters  & Processing
source("Functions/Dataprocess/Loading.R")

cbbPalette <- Loading_Color()
bac <- Loading_Table("Input/taxa_table 2.xls")
group <- Loading_Group("Input/group.txt")
group$Group <- factor(group$Group)
tree <- Loading_Tree("Input/phylogeny.tre")

result <- Data_Process(bac,group,tree)
bac1 <- result[[1]]
bac2 <- result[[2]]
bac3 <- result[[3]]
tree <- result[[4]]
Group_numb <- result[[5]]
Sample_numb <- result[[6]]
wid <- result[[7]]

## network analysis
dir.create("Results")
dir.create("Results/10.Network")
source("Functions/Network/network.calculation.R")
source("Functions/Network/cohesion.R")

### Group
cohes <- c()
for (i in 1:length(levels(group$Group))) {
    dir.create(paste("Results/10.Network/",levels(group$Group)[i],sep = ""))
    dir.create(paste("Results/10.Network/",levels(group$Group)[i],
                     "/01.Network_results",sep = ""))
    otu <- bac2[,group[group$Group == levels(group$Group)[i],"variable"]]
    otu2 <- bac3[,c(group[group$Group == levels(group$Group)[i],"variable"],
                    "Phylum","Class","Order","Family","Genus","Species")]
    result <- network.calculation(otu,otu2)
    write.table(result[[1]],
                paste("Results/10.Network/",levels(group$Group)[i],
                "/03.KeystoneTaxa/adjacency_unweight.txt",sep = ""),
                sep="\t",quote=F,col.names=NA)
    write.table(result[[2]],
                paste("Results/10.Network/",levels(group$Group)[i],
                "/01.Network_results/edge.csv",sep = ""),
                sep=",",quote=F,row.names = FALSE)
    coh <- cohesion.network(network_data1)
    
    neg <- coh$`Negative Cohesion`
    pos <- coh$`Positive Cohesion`
    neg.pos.coh <- data.frame(variable = rownames(pos),value = -neg/pos)
    cohes <- as.data.frame(rbind(cohes,neg.pos.coh))
}

dir.create("Results/10.Network/Stability_Group")
Group <- c()
for (i in 1: length(levels(group$Group))) {
    aa <- rep(levels(group$Group)[i],40)
    Group <- c(Group,aa)
}

### Cohesion
neg.pos.coh <- merge(cohes,group)
write.table(neg.pos.coh,"Results/10.Network/Stability_Group/cohesion.txt",
            sep = "\t",row.names = FALSE)
source("Functions/Network/cohesion.group.R")
result <- cohesion.group(neg.pos.coh,group)
pdf("Results/10.Network/Stability_Group/cohesion.pdf",
    width = 0.7 + Group_numb*0.5,height = 3.5)
result[[1]]
dev.off()
sink("Results/10.Network/Stability_Group/cohesion_test.txt")
result[[2]]
sink()

