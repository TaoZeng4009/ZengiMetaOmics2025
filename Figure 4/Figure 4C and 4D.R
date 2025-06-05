library(phyloseq)
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
pdf("Results/12.NullModel/Deter_Stoch.pdf",width = 0.7 + Group_numb*0.8,height = 4)
result[[5]]
dev.off()
write.table(result[[6]],"Results/12.NullModel/bNTI_RCbary.txt",
            sep = "\t",row.names = FALSE)
write.table(result[[7]],"Results/12.NullModel/sdratio.txt",
            sep = "\t",row.names = FALSE)

### Different of phylogenetic distance
dist <- result[[6]]
dist <- dist[dist$Group_1 == dist$Group_2,]
dist <- dist[,c(3,7)]
colnames(dist) <- c("Num","Group")
dist$Group <- factor(dist$Group)

x <- c("a","b")
y <- c("a","a")
if (length(levels(dist$Group)) == 2) {
    fit1 <- t.test(Num~Group,data = dist)
    dd <- dist %>%
        group_by(Group) %>%
        summarise(Max = max(Num))
    test <- data.frame(Group = levels(dist$Group),
                       value.x = if(fit1$p.value < 0.05){
                           x
                       }else{
                           y
                       },
                       value.y = dd$Max*1.03)
}else{
    fit1 <- aov(Num~Group,data = dist)
    tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
    res1 <- cld(tuk1,level = 0.05)
    dd <- dist %>%
        group_by(Group) %>%
        summarise(Max = max(Num))
    test <- data.frame(Group = levels(dist$Group),
                       value.x = res1$mcletters$Letters,
                       value.y = dd$Max*1.03)
}

p <- ggplot(dist,aes(Group,Num)) + 
    geom_boxplot(aes(color = Group),outlier.size = 0,size = 0.8) +
    geom_jitter(aes(color = Group),width = 0.25,alpha = 0.5) + 
    geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
              size = 8,color = "black") +
    scale_color_manual(values = cbbPalette) +
    ylab("Phylogenetic distance (betaMNTD)") +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(colour='black', size=18,
                                    vjust = 1.5),
          axis.text.y=element_text(colour='black',size=13),
          axis.text.x=element_text(colour = "black",size = 14),
          legend.position = "none")

pdf("Results/12.NullModel/Diff_bNMTD.pdf",width = 0.9 + 1.1*Group_numb,height = 5.5)
p
dev.off()
