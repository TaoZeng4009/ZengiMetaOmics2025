library(ggplot2)

dir.create("Results")
cbbPalette <- c("#B2182B","#56B4E9","#009E73")

## PCoA
library(vegan)
library(ape)
dir.create("Results/02.PCoA")
Gene <- read.table("Input/genes.fpkm.txt",header = TRUE,sep = "\t",row.names = 1)
group <- read.table("Input/group.txt",header = TRUE,sep = "\t")
colnames(group) <- c("sample","Group")
group$Group <- factor(group$Group,levels = c("In-TC","Ex-TC","QY"))
meta <- Gene
ARG_beta <- t(meta[,group[,1]])
ARG_dis <- vegdist(ARG_beta)

pcoa<- pcoa(ARG_dis, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2)
colnames(plotdata) <-c("sample","PC1","PC2")
plotdata <- merge(plotdata,group)
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)

ARG.adonis <- adonis2(ARG_dis~Group,data = plotdata)

ARG_pcoa2<-ggplot(plotdata, aes(PC1, PC2)) +
    geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
    stat_ellipse(aes(fill = Group),geom = "polygon",level = 0.95,alpha = 0.3)+
    geom_vline(aes(xintercept = 0),linetype="dotted")+
    geom_hline(aes(yintercept = 0),linetype="dotted")+
    geom_text(aes(0,-0.4),label = paste("Adonis test: R2 = ",round(ARG.adonis$R2[1],3),
                                         ", p-value = ",round(ARG.adonis$`Pr(>F)`[1],3),sep = ""),
              size = 5) + 
    scale_fill_manual(values=cbbPalette)+
    labs(title="PCoA - Lung transcriptome") + 
    xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
    ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
    theme(text=element_text(size=18))+
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(), 
          axis.title = element_text(color='black',size=18),
          axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_text(colour='black', size=18),
          axis.title.y=element_text(colour='black', size=18),
          axis.text=element_text(colour='black',size=16),
          legend.title=element_text(size = 12,face = "bold"),
          legend.text=element_text(size=12),
          legend.key=element_blank(),
          legend.background = element_rect(colour = "black"),
          legend.position = c(0.13,0.83))+
    theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))

pdf(file = "Results/02.PCoA/PCoA_ellipse.pdf",width = 5.4,height = 4.5)
ARG_pcoa2
dev.off()

#### Enrichment analysis
dir.create("Results/03.Enrichment")
enrichment <- read.table("Input/enrichment.txt",header = TRUE,sep = "\t")
enrichment$Rich_factor <- enrichment$Input.number/enrichment$Background.number
enrichment$Group <- factor(enrichment$Group,levels = c("In-TC-vs-Ex-TC",
                                                       "In-TC-vs-QY",
                                                       "Ex-TC-vs-QY"),
                           labels = c("In-TC\n-vs-\nEx-TC",
                                      "In-TC\n-vs-\nQY",
                                      "Ex-TC\n-vs-\nQY"))

p <-ggplot(enrichment,aes(Rich_factor,layer3)) +
    geom_point(aes(size = Input.number,color = Group)) +
    labs(title="Statistics of pathway enrichment for DEGs in lung",
         size = "DEGs") +
    scale_size_continuous(breaks = c(5,10)) +
    scale_color_manual(values = cbbPalette,guide = "none") + 
    facet_grid(Group~.,scales = "free_y",space = "free") + 
    xlab("Rich factor") + 
    ylab("") +
    theme(text=element_text(size=12))+
    theme(panel.background = element_rect(fill = "white",colour='black'),
          panel.grid.major = element_line(color = "grey",linetype = "dotted",linewidth = 0.3),
          panel.grid.minor = element_line(color = "grey",linetype = "dotted",linewidth = 0.3),
          axis.title = element_text(color='black',size=18),
          axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_text(colour='black', size=18),
          axis.text.x =element_text(colour='black',size=12),
          axis.text.y = element_text(colour='black',size=16),
          legend.title=element_text(size = 12,colour = "black"),
          legend.text=element_text(size=10),
          legend.key=element_blank(),legend.position = c(0.13,0.23),
          legend.background = element_rect(colour = "black"),
          legend.key.height=unit(0.4,"cm"),
          strip.text.y = element_text(angle = 0,colour = "black",size = 12)) +
    theme(plot.title = element_text(size=20,colour = "black",hjust = 1))

pdf("Results/03.Enrichment/enrichment.pdf",width = 10,height = 6)
p
dev.off()

### DAMs
dir.create("Results/04.Venn")
up <- read.table("Input/up.txt",header = TRUE,sep = "\t")
up <- list(`In-TC-vs-Ex-TC` = up$In.TC.vs.Ex.TC,
           `In-TC-vs-QY` = up$In.TC.vs.QY,
           `Ex-TC-vs-QY` = up$Ex.TC.vs.QY)
up$`Ex-TC-vs-QY` <- up$`Ex-TC-vs-QY`[up$`Ex-TC-vs-QY` != ""]
up$`In-TC-vs-Ex-TC` <- up$`In-TC-vs-Ex-TC`[up$`In-TC-vs-Ex-TC` != ""]
library(VennDiagram)
venn.diagram(up,filename = "Results/04.Venn/up.png",
             height = 5400,width = 5400,
             resolution = 600,imagetype = "png",units = "px",
             lwd = 2,lty = 1,fill = cbbPalette,cex = 1.5,
             cat.cex = 2,alpha = 0.8,margin = 0.05,fontface = 2,
             cat.fontface = 2, print.mode = c("raw","percent"))
inter <- get.venn.partitions(up)
aa <- unlist(inter$..values..[5])

down <- read.table("Input/down.txt",header = TRUE,sep = "\t")
down <- list(`In-TC-vs-Ex-TC` = down$In.TC.vs.Ex.TC,
           `In-TC-vs-QY` = down$In.TC.vs.QY,
           `Ex-TC-vs-QY` = down$Ex.TC.vs.QY)
down$`Ex-TC-vs-QY` <- down$`Ex-TC-vs-QY`[down$`Ex-TC-vs-QY` != ""]
down$`In-TC-vs-Ex-TC` <- down$`In-TC-vs-Ex-TC`[down$`In-TC-vs-Ex-TC` != ""]
venn.diagram(down,filename = "Results/04.Venn/down.png",
             height = 5400,width = 5400,
             resolution = 600,imagetype = "png",units = "px",
             lwd = 2,lty = 1,fill = cbbPalette,cex = 1.5,
             cat.cex = 2,alpha = 0.8,margin = 0.05,fontface = 2,
             cat.fontface = 2, print.mode = c("raw","percent"))
inter <- get.venn.partitions(down)
bb <- unlist(inter$..values..[5])

## Network
dir.create("Results/05.Network")
group <- read.table(file = "Input/group.txt",header = TRUE,sep = "\t")
otu <- read.table("Input/taxa_table.xls",header = FALSE,sep = "\t",quote = "",
                  row.names = 1)
colnames(otu) <- otu[1,]
otu <- otu[-1,]
for (i in 1:(ncol(otu)-1)) {
    otu[,i] <- as.numeric(otu[,i])
}
otu.taxonomy <- data.frame(V1 = rownames(otu),taxonomy = otu$taxonomy)
colnames(otu) <- gsub("-","_",colnames(otu))
otu <- otu[,c(group$variable,"taxonomy")]
otu <- otu[rowSums(otu[,1:(ncol(otu)-1)]) > 0,]

Gene <- Gene[rownames(Gene) %in% c(aa,bb),]
Gene <- Gene[,group$variable]

otu$taxonomy <- gsub(".*g__","",otu$taxonomy)
otu$taxonomy <- gsub(";.*","",otu$taxonomy)

library(linkET)
network_data <- cbind(t(otu[,-ncol(otu)]),t(Gene))
res <- fast_correlate(t(otu[,-ncol(otu)]),t(Gene),method = "spearman")

library(reshape2)
res1 <- res
data.r <- melt(res$r)
data.p <- melt(res$p)
result <- cbind(data.r,data.p)
result <- result[,c(1,2,3,6)]
colnames(result) <- c("Var1","Var2","R","P")

result.1 <- result[result$R > 0.8,]
result.2 <- result[result$R < -0.8,]
result <- rbind(result.1,result.2)
result <- result[result$P < 0.05,]
result<- result[result$R < 1,]
result$R[result$R > 0.8] = 1
result$R[result$R < -0.8] = -1
colnames(result) <- c("Source","Target","R","P")
edge <- result

nodes <- res %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    as_tibble(what = "vertices")
nodes$Type <- ifelse(nodes$name %in% aa,"Up-DEGs",
                     ifelse(nodes$name %in% bb,"Down-DEGs","Gut bacteria"))

library(tidyverse)
library(dplyr)
library(ggraph)
library(tidygraph)
net2 <- res %>%
    as_tbl_graph(abs(r) > 0.8, p < 0.05) %>%
    activate("nodes") %>%
    mutate(group = nodes$Type,
           Degree = centrality_degree())

p3 <- ggraph(net2, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Degree, shape = group),
                    show.legend = TRUE) +
    geom_node_text(aes(label = name),size = 1.5) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Degree") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/05.Network/Network_label.pdf",width = 10,height = 10)
p3
dev.off()

p4 <- ggraph(net2, "circular") +
    geom_edge_fan(color = "grey70", show.legend = FALSE, alpha = 0.2) +
    geom_node_point(aes(colour = group, size = Degree, shape = group),
                    show.legend = TRUE) +
    scale_edge_color_gradient2() +
    scale_size_continuous(name = "Degree)") +
    coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_blank())

pdf("Results/05.Network/Network_nolabel.pdf",width = 10,height = 10)
p4
dev.off()

## Correlation
dir.create("Results/06.Correlation")
aa <- unique(result$Source)
otu <- otu[aa,]
otu$taxonomy <- gsub(" .*","",otu$taxonomy)
rownames(otu) <- gsub("OTU","ASV",rownames(otu))
rownames(otu) <- paste(rownames(otu),"(",otu$taxonomy,")",sep = "")
bb <- c("CCL1","CCL19","CCR2","CCR4","CCR5","CCR6","CCR7","CCR9","CCR10","CSF1R",
        "CXCL13","CXCL13L2","CXCL13L3","CXCR5","CXCR7","IL18","IL18R1","IL2RA",
        "IL6","IL8L1","IL8L2","LOC100857191","TNFRSF10B","TNFRSF1A","XCR1")
Gene <- Gene[rownames(Gene) %in% bb,]

otu <- otu[,-ncol(otu)]
library(pheatmap)
library(psych)
res <- corr.test(t(otu),t(Gene),method = "spearman",adjust = "holm")
pdf("Results/06.Correlation/heatmap.pdf",width = 10,height = 12)
pheatmap(res$r, fontsize_number=13,fontsize = 12,number_color = "white",
         display_numbers = matrix(ifelse(res$p <= 0.01, "++", 
                                         ifelse(res$p <= 0.05 ,"+"," ")), nrow(res$p)))
dev.off()
