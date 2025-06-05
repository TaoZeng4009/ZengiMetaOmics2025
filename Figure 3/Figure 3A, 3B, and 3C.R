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
library(VennDiagram)
library(RColorBrewer)
library(gplots)
library(plotrix)
library(ggalluvial)
library(d3Network)
library(plotly)

### Loading Parameters  & Processing
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

## Core analysis
dir.create("Results")
dir.create("Results/04.Core")

### Venn analysis by groups
dir.create("Results/04.Core/01.Venn_group")

#### Venn diagram
source("Functions/Core/Flower.g.R")
pdf(file = "Results/04.Core/01.Venn_group/Venn_group.pdf",
    width = 8,height = 8)
flower.g(bac2,group)
dev.off()

source("Functions/Core/Core.species.id.g.R")
core_species_id_g <- core.species.id.g(bac2,group)
core_species_id_g <- bac3[core_species_id_g,(ncol(bac3)-5):ncol(bac3)]
core_species_id_g$ID <- rownames(core_species_id_g)
core_species_id_g <- core_species_id_g[,c("ID",colnames(core_species_id_g)[1:6])]
write.table(core_species_id_g,"Results/04.Core/01.Venn_group/Shared_species_group.txt",
            sep = "\t",row.names = FALSE)

#### Differences in total abundance of shared species among different groups
dir.create("Results/04.Core/02.Diff_group")
source("Functions/Core/Core.species.abun.g.R")
result <- core.species.abun.g(bac2,core_species_id_g,group)
pdf(file = "Results/04.Core/02.Diff_group/Shared_species_abundance_diff_jitter.pdf",
    width = wid*Group_numb,height = 4.2)
result[[2]]
dev.off()
sink("Results/04.Core/02.Diff_group/Shared_species_abundance_diff_test.txt")
summary(result[[4]])
sink()

#### Ratio of different phyla in shared species
dir.create("Results/04.Core/03.TAX_group")

### Sankey of shared species
source("Functions/Core/Sankey.core.species.R")
plotdata <- sankey.core.species(bac3,core_species_id_g)
nodes <- data.frame(name = unique(c(as.character(plotdata$source),
                                    as.character(plotdata$target))),
                    stringsAsFactors = FALSE)
nodes$ID <- 0:(nrow(nodes)-1)
sankey <- merge(plotdata,nodes,by.x = "source",by.y = "name")
sankey <- merge(sankey,nodes,by.x = "target",by.y = "name")
colnames(sankey) <- c("X","Y","value","source","target")
sankey <- subset(sankey,select = c("source","target","value"))
nodes <- subset(nodes,select = c("name"))

d3Sankey(Links = sankey,Nodes = nodes,fontsize = 20,
         Source = "source",Target = "target",Value = "value",
         NodeID = "name",fontFamily = "Arial",
         file = "Results/04.Core/03.TAX_group/Sankey.html",
         width = 1200,height = 900)
write.table(plotdata,"Results/04.Core/03.TAX_group/sankey_data.txt",
            row.names = FALSE,sep = "\t")
