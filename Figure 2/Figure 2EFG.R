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
### Betapart
dir.create("Results/14.Niche/04.betapart")

#### Unweighted
bac.u <- t(bac2)
bac.u[bac.u > 0] = 1
beta.bac.u <- betapart.core(bac.u)
bac.u.jac <- beta.pair(beta.bac.u,index.family = "jac")
bac.u.turn <- as.matrix(bac.u.jac$beta.jtu)
bac.u.turn[lower.tri(bac.u.turn)] = 0
bac.u.turn <- melt(bac.u.turn)
bac.u.turn <- bac.u.turn[bac.u.turn$value > 0,]

bac.u.nest <- as.matrix(bac.u.jac$beta.jne)
bac.u.nest[lower.tri(bac.u.nest)] = 0
bac.u.nest <- melt(bac.u.nest)
bac.u.nest <- bac.u.nest[bac.u.nest$value > 0,]

bac.u <- rbind(bac.u.turn,bac.u.nest)
bac.u$Type <- c(rep("Turnover",nrow(bac.u.turn)),
                rep("Nestedness",nrow(bac.u.nest)))
write.table(bac.u,"Results/14.Niche/04.betapart/betapart_Un_all.txt",
            sep = "\t",row.names = FALSE)

if (nrow(bac.u.nest) > 0) {
    x <- c("a","b")
    y <- c("a","a")
    bac.u$Type <- factor(bac.u$Type)
    test <- c()
    fit1 <- t.test(value~Type,data = bac.u)
    if (fit1$p.value < 0.05) {
        test <- x
    }else{
        test <- y
    }
    test <- data.frame(Type = levels(bac.u$Type),value = test)
    test.beta.u1 <- bac.u %>% 
        group_by(Type) %>% 
        summarise(Max = max(value))
    test <- merge(test,test.beta.u1)
}else{
    test <- data.frame(Type = c("Nestedness","Turnover"),
                       value = c("ND",""),
                       Max = c(0,0))
}

beta.u1 <- ggplot(bac.u,aes(Type,value)) + 
    geom_boxplot(aes(color = Type),outlier.size = 0,size = 0.8,width = .6) +
    geom_jitter(aes(color = Type),size = 1,alpha = 0.5,
                position = position_jitter(width = 0.3)) +
    geom_text(data = test,aes(x = Type,y = Max + max(Max)*0.1,label = value),
              size = 6,color = "black") +
    scale_color_manual(values = cbbPalette) +
    ylab("Partitioning beta diversity of species richness") +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(colour='black', size=11,
                                    vjust = 2,hjust = 1),
          axis.text.y=element_text(colour='black',size=10),
          axis.text.x=element_text(colour = "black",size = 12,
                                   angle = 45,hjust = 1,vjust = 1),
          legend.position = "none")

pdf("Results/14.Niche/04.betapart/Un_betapart_all.pdf",
    width = 2,height = 4)
beta.u1
dev.off()

#### Weighted
bac.w <- t(bac2)
beta.bac.w <- betapart.core.abund(bac.w)
bac.w.jac <- beta.pair.abund(beta.bac.w)
bac.w.turn <- as.matrix(bac.w.jac$beta.bray.bal)
bac.w.turn[lower.tri(bac.w.turn)] = 0
bac.w.turn <- melt(bac.w.turn)
bac.w.turn <- bac.w.turn[bac.w.turn$value > 0,]

bac.w.nest <- as.matrix(bac.w.jac$beta.bray.gra)
bac.w.nest[lower.tri(bac.w.nest)] = 0
bac.w.nest <- melt(bac.w.nest)
bac.w.nest <- bac.w.nest[bac.w.nest$value > 0,]

bac.w <- rbind(bac.w.turn,bac.w.nest)
bac.w$Type <- c(rep("Turnover",nrow(bac.w.turn)),
                rep("Nestedness",nrow(bac.w.nest)))
write.table(bac.u,"Results/14.Niche/04.betapart/betapart_W_all.txt",
            sep = "\t",row.names = FALSE)

if (nrow(bac.w.nest) > 0) {
    x <- c("a","b")
    y <- c("a","a")
    bac.w$Type <- factor(bac.w$Type)
    test <- c()
    fit1 <- t.test(value~Type,data = bac.w)
    if (fit1$p.value < 0.05) {
        test <- x
    }else{
        test <- y
    }
    test <- data.frame(Type = levels(bac.w$Type),value = test)
    test.beta.w1 <- bac.w %>% 
        group_by(Type) %>% 
        summarise(Max = max(value))
    test <- merge(test,test.beta.w1)
}else{
    test <- data.frame(Type = c("Nestedness","Turnover"),
                       value = c("ND",""),
                       Max = c(0,0))
}


beta.w1 <- ggplot(bac.w,aes(Type,value)) + 
    geom_boxplot(aes(color = Type),outlier.size = 0,size = 0.8,width = .6) +
    geom_jitter(aes(color = Type),size = 1,alpha = 0.5,
                position = position_jitter(width = 0.3)) +
    geom_text(data = test,aes(x = Type,y = Max + max(Max)*0.1,label = value),
              size = 5,color = "black") +
    scale_color_manual(values = cbbPalette) +
    ylab("Partitioning beta diversity of species abundances") +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(colour='black', size=11,
                                    vjust = 2,hjust = 1),
          axis.text.y=element_text(colour='black',size=10),
          axis.text.x=element_text(colour = "black",size = 12,
                                   angle = 45,hjust = 1,vjust = 1),
          legend.position = "none")

pdf("Results/14.Niche/04.betapart/W_betapart_all.pdf",
    width = 2,height = 4)
beta.w1
dev.off()

### Rao quadratic entropy
dir.create("Results/14.Niche/06.Rao")
source("Functions/Niche/f(rao).R")
tree <- compute.brlen(tree, power = 10)
bb <- c()
for (i in 1:length(levels(group$Group))) {
    bac.g <- bac2[,group$variable[group$Group == levels(group$Group)[i]]]
    tree.g <- prune.sample(t(bac.g),tree)
    bac.g <- bac.g[tree.g$tip.label,]
    result <- Rao(community.composition = bac.g,phylo.distance = tree.g,
                  funct.distance = NULL,local.community.ab.weight = FALSE,
                  outputlike = "tables")
    aa <- data.frame(Group = rep(levels(group$Group)[i],
                                 nrow(result$PD$Pairwise.partition.Phylogeny)),
                     Q.Beta.st = result$PD$Pairwise.partition.Phylogeny[,"Q.Beta.st"],
                     Q.alpha.st = 1-result$PD$Pairwise.partition.Phylogeny[,"Q.Beta.st"])
    bb <- as.data.frame(rbind(bb,aa))
}

write.table(bb,"Results/14.Niche/06.Rao/Rao.txt",sep = "\t",row.names = FALSE)

bb <- melt(bb)

bb$Group <- factor(bb$Group,levels = levels(group$Group))
alpha.test <- bb[bb$variable == "Q.alpha.st",]
if (length(levels(group$Group)) == 2) {
    fit1 <- t.test(value~Group,data = alpha.test)
    dd <- alpha.test %>%
        group_by(Group) %>%
        summarise(Max = max(value))
    test <- data.frame(Group = levels(alpha.test$Group),
                       value.x = if(fit1$p.value < 0.05){
                           x
                       }else{
                           y
                       },
                       value.y = dd$Max*1.03)
}else{
    fit1 <- aov(value~Group,data = alpha.test)
    tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
    res1 <- cld(tuk1,level = 0.05)
    dd <- alpha.test %>%
        group_by(Group) %>%
        summarise(Max = max(value))
    test <- data.frame(Group = levels(alpha.test$Group),
                       value.x = res1$mcletters$Letters,
                       value.y = dd$Max*1.03)
    test$variable <- rep("Alpha",nrow(test))
}

aa <- c()
if (length(levels(alpha.test$Group)) == 2) {
    aa <- c(aa,fit1$p.value)
    test.result <- data.frame(Index = colnames(alpha1)[2:ncol(alpha1)],
                              pvalue = aa)
}else{
    aa <- rbind(aa,summary(tuk1)$test$pvalue)
    colnames(aa) <- names(summary(tuk1)$test$coefficients)
    test.result <- as.data.frame(t(aa))
    colnames(test.result) <- "Q.alpha.st"
    test.result$Index <- rownames(test.result)
}

alpha.test <-  bb[bb$variable == "Q.Beta.st",]
if (length(levels(group$Group)) == 2) {
    fit1 <- t.test(value~Group,data = alpha.test)
    dd <- alpha.test %>%
        group_by(Group) %>%
        summarise(Max = max(value))
    test2 <- data.frame(Group = levels(alpha.test$Group),
                        value.x = if(fit1$p.value < 0.05){
                            x
                        }else{
                            y
                        },
                        value.y = dd$Max*1.03)
}else{
    fit1 <- aov(value~Group,data = alpha.test)
    tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
    res1 <- cld(tuk1,level = 0.05)
    dd <- alpha.test %>%
        group_by(Group) %>%
        summarise(Max = max(value))
    test2 <- data.frame(Group = levels(alpha.test$Group),
                        value.x = res1$mcletters$Letters,
                        value.y = dd$Max*1.03)
    test2$variable <- rep("Beta",nrow(test2))
}

aa <- c()
if (length(levels(alpha.test$Group)) == 2) {
    aa <- c(aa,fit1$p.value)
    test.result <- data.frame(Index = colnames(alpha1)[2:ncol(alpha1)],
                              pvalue = aa)
}else{
    aa <- rbind(aa,summary(tuk1)$test$pvalue)
    colnames(aa) <- names(summary(tuk1)$test$coefficients)
    test.result2 <- as.data.frame(t(aa))
    colnames(test.result2) <- "Q.beta.st"
    test.result2$Index <- rownames(test.result2)
}
test.result <- merge(test.result,test.result2)
    

write.table(test.result2,"Results/14.Niche/06.Rao/Rao_test.txt",
            sep = "\t",row.names = FALSE)

test <- as.data.frame(rbind(test,test2))

dd <- max(str_length(levels(group$Group)))

bb$variable <- gsub("Q.Beta.st","Beta",bb$variable)
bb$variable <- gsub("Q.alpha.st","Alpha",bb$variable)
bb$value <- bb$value*100
test$value.y <- test$value.y*100
beta.u1 <- ggplot(bb,aes(Group,value)) + 
    geom_boxplot(aes(color = Group),outlier.size = 0,size = 0.8,width = .6) +
    geom_jitter(aes(color = Group),size = 1,alpha = 0.5,
                position = position_jitter(width = 0.3)) +
    geom_text(data = test,aes(x = Group,y = value.y + max(value.y)*0.001,label = value.x),
              size = 5,color = "black") +
    scale_color_manual(values = cbbPalette) +
    ylab("Contribution (%) to gamma diversity") +
    facet_grid(variable~.,scales = "free_y") +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(colour='black', size=11),
          axis.text.y=element_text(colour='black',size=10),
          axis.text.x=element_text(colour = "black",size = 12,
                                   angle = ifelse(dd > 4,45,0),
                                   hjust = ifelse(dd > 4,1,0.5),
                                   vjust = 1),
          strip.text = element_text(colour = "black",size = 14),
          legend.position = "none")

pdf("Results/14.Niche/06.Rao/Rao.pdf",
    width = 1 + Group_numb*0.5,height = 3.5)
beta.u1
dev.off()
