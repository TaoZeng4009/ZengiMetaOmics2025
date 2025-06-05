library(ggplot2)
library(multcomp)
library(tidyverse)
library(reshape2)
cbbPalette <- c("#B2182B","#56B4E9","#E69F00","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999",
                "#ADD1E5")
data <- read.table("data.txt",header = TRUE,sep = "\t")
data$Group <- factor(data$Group,levels = c("In-TC","Ex-TC","BY","CH","LS","QY"))

alpha.phy <- data[,1:4]
aa <- c()
for (i in 2:ncol(alpha.phy)) {
    alpha.test <- alpha.phy[,c(i,1)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max + max(dd$Max)*0.05,
                           variable = rep(colnames(alpha.phy)[i],2))
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max + max(dd$Max)*0.05,
                           variable = rep(colnames(alpha.phy)[i],length(levels(alpha.test$Group))))
    }
    aa <- as.data.frame(rbind(aa,test))
}

test <- aa
alpha.test <- alpha.phy
alpha.test <- melt(alpha.test)

alpha.test$variable <- factor(alpha.test$variable)
alpha.test$Group <- factor(alpha.test$Group,levels = levels(data$Group))
test$Group <- factor(test$Group,levels = levels(data$Group))

p <- ggplot(alpha.test,aes(Group,value,color = Group)) + 
    geom_boxplot(width = 0.6,outlier.color = "transparent") +
    geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
    geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
              size = 4.5,color = "black") +
    labs(y = "Enzyme activities",
         x = "") +
    scale_color_manual(values = cbbPalette) +
    facet_wrap(variable~.,scales = "free_y",nrow = 1) + 
    theme_bw()+
    theme(panel.grid=element_blank(),
          axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y = element_text(color = "black",size = 16),
          axis.text.y=element_text(colour='black',size=10),
          axis.text.x=element_text(colour = "black",size = 14,
                                   angle = 45,vjust = 1,hjust = 1),
          legend.position = "none",
          strip.text = element_text(colour = "black",size = 14),
          strip.text.y = element_text(colour = "black",size = 14))

pdf("Enzyme.activities.pdf",height = 2.5,width = 7)
p
dev.off()

alpha.phy <- data[,c(1,5:8)]
aa <- c()
for (i in 2:ncol(alpha.phy)) {
    alpha.test <- alpha.phy[,c(i,1)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max + max(dd$Max)*0.05,
                           variable = rep(colnames(alpha.phy)[i],2))
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max + max(dd$Max)*0.05,
                           variable = rep(colnames(alpha.phy)[i],length(levels(alpha.test$Group))))
    }
    aa <- as.data.frame(rbind(aa,test))
}

test <- aa
alpha.test <- alpha.phy
alpha.test <- melt(alpha.test)

alpha.test$variable <- factor(alpha.test$variable)
alpha.test$Group <- factor(alpha.test$Group,levels = levels(data$Group))
test$Group <- factor(test$Group,levels = levels(data$Group))

p <- ggplot(alpha.test,aes(Group,value,color = Group)) + 
    geom_boxplot(width = 0.6,outlier.color = "transparent") +
    geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
    geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
              size = 4.5,color = "black") +
    labs(y = "Myocardium parameters",
         x = "") +
    scale_color_manual(values = cbbPalette) +
    facet_wrap(variable~.,scales = "free_y",nrow = 2) + 
    theme_bw()+
    theme(panel.grid=element_blank(),
          axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y = element_text(color = "black",size = 20),
          axis.text.y=element_text(colour='black',size=10),
          axis.text.x=element_text(colour = "black",size = 14,
                                   angle = 45,vjust = 1,hjust = 1),
          legend.position = "none",
          strip.text = element_text(colour = "black",size = 14),
          strip.text.y = element_text(colour = "black",size = 14))

pdf("Myocardium.pdf",height = 4.5,width = 5)
p
dev.off()

alpha.phy <- data[,c(1,9:14)]
aa <- c()
for (i in 2:ncol(alpha.phy)) {
    alpha.test <- alpha.phy[,c(i,1)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max + max(dd$Max)*0.05,
                           variable = rep(colnames(alpha.phy)[i],2))
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max + max(dd$Max)*0.05,
                           variable = rep(colnames(alpha.phy)[i],length(levels(alpha.test$Group))))
    }
    aa <- as.data.frame(rbind(aa,test))
}

test <- aa
alpha.test <- alpha.phy
alpha.test <- melt(alpha.test)

alpha.test$variable <- factor(alpha.test$variable)
alpha.test$Group <- factor(alpha.test$Group,levels = levels(data$Group))
test$Group <- factor(test$Group,levels = levels(data$Group))

p <- ggplot(alpha.test,aes(Group,value,color = Group)) + 
    geom_boxplot(width = 0.6,outlier.color = "transparent") +
    geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
    geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
              size = 4.5,color = "black") +
    labs(y = "Indices of hematology",
         x = "") +
    scale_color_manual(values = cbbPalette) +
    facet_wrap(variable~.,scales = "free_y",nrow = 2) + 
    theme_bw()+
    theme(panel.grid=element_blank(),
          axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y = element_text(color = "black",size = 20),
          axis.text.y=element_text(colour='black',size=10),
          axis.text.x=element_text(colour = "black",size = 14,
                                   angle = 45,vjust = 1,hjust = 1),
          legend.position = "none",
          strip.text = element_text(colour = "black",size = 14),
          strip.text.y = element_text(colour = "black",size = 14))

pdf("Hematology.pdf",height = 4.5,width = 7)
p
dev.off()

alpha.phy <- data[,c(1,15:17)]
aa <- c()
for (i in 2:ncol(alpha.phy)) {
    alpha.test <- alpha.phy[,c(i,1)]
    colnames(alpha.test) <- c("Num","Group")
    alpha.test$Group <- factor(alpha.test$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(alpha.test$Group)) == 2) {
        fit1 <- t.test(Num~Group,data = alpha.test)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max + max(dd$Max)*0.05,
                           variable = rep(colnames(alpha.phy)[i],2))
    }else{
        fit1 <- aov(Num~Group,data = alpha.test)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1)
        dd <- alpha.test %>%
            group_by(Group) %>%
            summarise(Max = max(Num))
        test <- data.frame(Group = levels(alpha.test$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max + max(dd$Max)*0.05,
                           variable = rep(colnames(alpha.phy)[i],length(levels(alpha.test$Group))))
    }
    aa <- as.data.frame(rbind(aa,test))
}

test <- aa
alpha.test <- alpha.phy
alpha.test <- melt(alpha.test)

alpha.test$variable <- factor(alpha.test$variable)
alpha.test$Group <- factor(alpha.test$Group,levels = levels(data$Group))
test$Group <- factor(test$Group,levels = levels(data$Group))

p <- ggplot(alpha.test,aes(Group,value,color = Group)) + 
    geom_boxplot(width = 0.6,outlier.color = "transparent") +
    geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
    geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
              size = 4.5,color = "black") +
    labs(y = "Weight indices",
         x = "") +
    scale_color_manual(values = cbbPalette) +
    facet_wrap(variable~.,scales = "free_y",nrow = 1) + 
    theme_bw()+
    theme(panel.grid=element_blank(),
          axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y = element_text(color = "black",size = 16),
          axis.text.y=element_text(colour='black',size=10),
          axis.text.x=element_text(colour = "black",size = 14,
                                   angle = 45,vjust = 1,hjust = 1),
          legend.position = "none",
          strip.text = element_text(colour = "black",size = 14),
          strip.text.y = element_text(colour = "black",size = 14))

pdf("Weight.pdf",height = 2.5,width = 7)
p
dev.off()