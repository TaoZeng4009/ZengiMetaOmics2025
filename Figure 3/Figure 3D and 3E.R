ibrary(phyloseq)
library(tidyverse)
library(vegan)
library(picante)
library(minpack.lm)
library(FSA)
library(eulerr)
library(ggplot2)
library(grid)
require(Hmisc)
require(stats4)
library(parallel)
library(ggstatsplot)
library(RColorBrewer)
library(reshape2)
library(multcomp)

dir.create("Results")
dir.create("Results/12.NullModel")

cbbPalette <- c("#B2182B","#56B4E9","#E69F00","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",
                "#ADD1E5")

map = read.delim("Input/group.txt")
otu = read.delim("Input/taxa_table.xls", row.names=1,header = FALSE)
tree <- read.tree("Input/phylogeny.tre")

Group_numb <- length(unique(map$Group))

colnames(otu) <- otu[1,]
otu <- otu[-1,]
for (i in 1:(ncol(otu)-1)) {
    otu[,i] <- as.numeric(otu[,i])
}

rownames(map) <- map[,1]
map$Group <- factor(map$Group,levels = unique(map$Group))
otu <- otu[,c(map$variable,"taxonomy")]

tax <- otu %>% 
    separate(taxonomy,c("Domain","Phylum","Class","Order","Family","Genus","Species"),";")
tax <- tax[,(ncol(tax)-6):ncol(tax)]

otu <- otu[,-ncol(otu)]

source("Functions/NullModel/inputMicro.R")
ps = inputMicro(otu,tax,map,tree,group  = "Group")

source("Functions/NullModel/betaNTI.R")
result = bNTICul(ps = ps ,group  = "Group",num = 499,thread = 24)
bNTI = result[[1]]
write.csv(bNTI, "Results/12.NullModel/bNTI.csv")

source("Functions/NullModel/RCbray.R")
result = RCbary(ps = ps ,group  = "Group",num = 499,thread = 24)
RCbary = result[[1]]
write.csv(RCbary, "Results/12.NullModel/RCbary.csv")


source("Functions/NullModel/betaNTIplot.R")
bNTI = read.csv("Results/12.NullModel/bNTI.csv",row.names = 1)
RCb = read.csv("Results/12.NullModel/RCbary.csv",row.names = 1)

result = bNTIRCPlot(ps = ps ,RCb  = RCb,bNTI = bNTI,group  = "Group")

pdf("Results/12.NullModel/bNTI.pdf",width = 0.7 + Group_numb*0.8,height = 3)
result[[1]]
dev.off()
pdf("Results/12.NullModel/process.pdf",width = Group_numb*0.4 + 3.4,height = 3)
result[[2]]
dev.off()
