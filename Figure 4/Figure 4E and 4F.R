### Pipeline for standard analyses of alpha diversity of bacterial communities
### Based on R v4.0.2
library(vegan)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ape)
library(picante)
library(dplyr)
library(TITAN2)
library(linkET)
library(betapart)
library(multcomp)
library(SoDA)
library(phyloseq)
library(spaa)

## Parameters
source("Functions/Dataprocess/Loading.R")

cbbPalette <- Loading_Color()
bac <- Loading_Table("Input/taxa_table.xls")
group <- Loading_Group("Input/group.txt")
group$Group <- factor(group$Group,levels = unique(group$Group))
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
dir.create("Results/14.Niche")

### Functional redundancy
dir.create("Results/14.Niche/03.FunRedundancy")
funcr <- read.table("Input/Tax4Fun2_FRI/relative_functional_redundancy.txt",
                    sep = "\t",row.names = 1,header = FALSE)
colnames(funcr) <- funcr[1,]
funcr <- funcr[-1,]
for (i in 1:(ncol(funcr)-1)) {
    funcr[,i] <- as.numeric(funcr[,i])
}
funcr <- funcr[,group$variable]

site_dis <- as.matrix(vegdist(t(funcr),method = "euclidean"))
diag(site_dis) <- 0
site_dis[upper.tri(site_dis)] <- 0
site_dis <- melt(site_dis)
site_dis <- subset(site_dis, value != 0)

bMNTD <- read.table("Results/12.NullModel/bNTI_RCbary.txt",sep = "\t",header = TRUE)

dist <- bMNTD
dist <- dist[dist$Group_1 == dist$Group_2,]
data.p <- dist[,c(1,2,6,7)]
colnames(data.p) <- c("Var1","Var2","bNMTD","Group")

comm_dis.p <- merge(data.p,site_dis,by = c("Var1","Var2"))
colnames(comm_dis.p) <- c("site1","site2","x","Group","y")

p <- ggplot(comm_dis.p,aes(x,y,Group = Group)) +
    geom_point(aes(color = Group),size = 2,alpha = 0.7) + 
    geom_smooth(aes(color = Group,fill = Group),
                method = 'lm', formula = y ~ x) +
    scale_fill_manual(values = cbbPalette) + 
    scale_color_manual(values = cbbPalette) +
    theme_bw() + 
    labs(y = "ΔFunctional redundancy index",
         x = "βNTI") +
    theme(axis.title.x = element_text(size = 12,colour = "black"),
          axis.title.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 10,colour = "black"),
          axis.text.y = element_text(size = 10,colour = "black"),
          legend.title = element_text(size = 12,colour = "black"),
          legend.text = element_text(size = 10,colour = "black"))

pdf("Results/14.Niche/03.FunRedundancy/FunR_bNTI.pdf",width = 5.5,height = 4)
p
dev.off()
png("Results/14.Niche/03.FunRedundancy/FunR_bNTI.png",
    width = 3000,height = 2200,res = 600)
p
dev.off()

test <- c()
for (i in 1:length(unique(comm_dis.p$Group))) {
    m <- lm(y~x,comm_dis.p[comm_dis.p$Group == unique(comm_dis.p$Group)[i],])
    a <- summary(m)
    aa <- c(m$coefficients[2],m$coefficients[1],a$adj.r.squared,a$coefficients[2,4])
    test <- rbind(test,aa)
}
test <- as.data.frame(test)
colnames(test) <- c("slope","intercept","r2","pvalue")
rownames(test) <- unique(comm_dis.p$Group)
write.table(test,"Results/14.Niche/03.FunRedundancy/FunR_bNTI.txt",sep = "\t")

funr1 <- data.frame(Group1 = rowMeans(funcr[,group[group$Group == levels(group$Group)[1],"variable"]]),
                    Group2 = rowMeans(funcr[,group[group$Group == levels(group$Group)[2],"variable"]]))
funr1 <- funr1[funr1$Group1 + funr1$Group2 > 0,]
funr1$ratio <- log(funr1$Group1/funr1$Group2)
funr1$X <- 1:nrow(funr1)
funr1$ratio <- ifelse(funr1$ratio == Inf,10,
                      ifelse(funr1$ratio == -Inf,-10,funr1$ratio))
funr1$Group <- ifelse(funr1$ratio > 0,
                      paste(levels(group$Group)[1]," (n=",nrow(funr1[funr1$ratio > 0,]),")",sep = ""),
                      paste(levels(group$Group)[2]," (n=",nrow(funr1[funr1$ratio < 0,]),")",sep = ""))

p <- ggplot(funr1,aes(X,ratio)) + 
    geom_point(aes(color = Group)) + 
    scale_color_manual(values = cbbPalette) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(limits = c(min(funr1$ratio)*1.1,max(funr1$ratio)*1.1)) + 
    theme_bw() + 
    labs(x = "",
         y = paste("log ratio(FRI ",levels(group$Group)[1],"/FRI ",levels(group$Group)[2],")",sep = "")) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10,colour = "black"),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 10,colour = "black"))

pdf("Results/14.Niche/03.FunRedundancy/FRI.pdf",width = 5,height = 4)
p
dev.off()    
