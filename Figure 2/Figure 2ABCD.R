### Pipeline for standard analyses of beta diversity of bacterial communities
### Based on R v4.0.2
library(ggplot2)
library(reshape2)
library(vegan)
library(dplyr)
library(tidyverse)
library(tidyr)
library(patchwork)
library(multcomp)
library(RColorBrewer)
library(pheatmap)
library(gplots)
library(plotrix)
library(ape)
library(ggrepel)
library(FactoMineR)
library(picante)
library(scatterplot3d)
library(phyloseq)
library(GUniFrac)
library(corrplot)
library(grid)
library(patchwork)

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
species.type <- "Microbial"

### Beta diversity analysis
dir.create("Results/05.Beta")
dir.create("Results/05.Beta/04.PCoA")
dir.create("Results/05.Beta/08.Beta_test")

### Unifrac distance
unifrac <- phyloseq(otu_table(t(bac2),taxa_are_rows = F),phy_tree(tree))
W.unifrac <- distance(unifrac,method = "wunifrac")
U.unifrac <- distance(unifrac,method = "unifrac")
W.unifrac <- as.matrix(W.unifrac)
write.table(W.unifrac,
            "Results/05.Beta/01.beta_diversity/Weighted_Unifrac_distance.txt",
            sep = "\t")
U.unifrac <- as.matrix(U.unifrac)
write.table(U.unifrac,
            "Results/05.Beta/01.beta_diversity/Unweighted_Unifrac_distance.txt",
            sep = "\t")

#### PCoA of Weighted_Unifrac distance
dir.create("Results/05.Beta/04.PCoA/Weighted_Unifrac")
source("Functions/Beta/PCoA2.R")
result <- pcoa.community2(W.unifrac,group,species.type)
write.table(result[[1]],"Results/05.Beta/04.PCoA/Weighted_Unifrac/pcoa_top2.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/05.Beta/04.PCoA/Weighted_Unifrac/Pcoa_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
write.table(result[[5]],
            "Results/05.Beta/04.PCoA/Weighted_Unifrac/PCoA_score.txt",
            sep = "\t")

#### PCoA of Weighted_Unifrac distance with statistical tests
pdf(file = "Results/05.Beta/04.PCoA/Weighted_Unifrac/PCoA_with_tests.pdf",
    width = 14,height = 12)
result[[7]]
dev.off()

#### PCoA of Unweighted_Unifrac distance
dir.create("Results/05.Beta/04.PCoA/Unweighted_Unifrac")
source("Functions/Beta/PCoA2.R")
result <- pcoa.community2(U.unifrac,group,species.type)
write.table(result[[1]],"Results/05.Beta/04.PCoA/Unweighted_Unifrac/pcoa_top2.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/05.Beta/04.PCoA/Unweighted_Unifrac/Pcoa_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
write.table(result[[5]],
            "Results/05.Beta/04.PCoA/Unweighted_Unifrac/PCoA_score.txt",
            sep = "\t")

#### Distance heatmap for Weighted_Unifrac distance
dist_dd <- as.matrix(W.unifrac)

result <- diff.beta(dist_dd,group)
pdf(file = "Results/05.Beta/02.Distance/Weighted_Unifrac/Group_distance_diff_jitter.pdf",
    width = Group_numb*0.5 + 0.5,height = 4.2)
result[[2]]
dev.off()
sink("Results/05.Beta/02.Distance/Weighted_Unifrac/Group_distance_diff_test.txt")
summary(result[[4]])
sink()

#### Distance heatmap for Unweighted_Unifrac distance
dist_dd <- as.matrix(U.unifrac)

result <- diff.beta(dist_dd,group)
pdf(file = "Results/05.Beta/02.Distance/Unweighted_Unifrac/Group_distance_diff_jitter.pdf",
    width = Group_numb*0.5 + 0.5,height = 4.2)
result[[2]]
dev.off()
sink("Results/05.Beta/02.Distance/Unweighted_Unifrac/Group_distance_diff_test.txt")
summary(result[[4]])
sink()

### ANOSIM, Adonis, and MRPP tests of bary-curtis distance
dir.create("Results/05.Beta/08.Beta_test/Weighted_Unifrac")
result <- beta.test(W.unifrac,group)
sink("Results/05.Beta/08.Beta_test/Weighted_Unifrac/Adonis.txt")
result[[1]]
sink()

### ANOSIM, Adonis, and MRPP tests of bary-curtis distance
dir.create("Results/05.Beta/08.Beta_test/Unweighted_Unifrac")
result <- beta.test(U.unifrac,group)
sink("Results/05.Beta/08.Beta_test/Unweighted_Unifrac/Adonis.txt")
result[[1]]
sink()

