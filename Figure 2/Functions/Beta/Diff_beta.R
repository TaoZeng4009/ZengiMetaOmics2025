diff.beta <- function(ARG_sub,group){
    ARG_sub1 <- melt(ARG_sub)
    ARG_sub1 <- ARG_sub1[ARG_sub1$Var1 != ARG_sub1$Var2,]
    colnames(group) <- c("Var1","Group1")
    ARG_sub1 <- merge(ARG_sub1,group)
    colnames(group) <- c("Var2","Group2")
    ARG_sub1 <- merge(ARG_sub1,group)
    ARG_sub1 <- ARG_sub1[ARG_sub1$Group1 == ARG_sub1$Group2,]
    ARG_sub1$Group <- ARG_sub1$Group1
    core_ARG_ta <- ARG_sub1[,c(3,6)]
    colnames(core_ARG_ta) <- c("Abun","Group")
    core_ARG_ta$Group <- factor(core_ARG_ta$Group)
    
    if (length(levels(core_ARG_ta$Group)) == 2) {
        fit1 <- t.test(Abun~Group,data = core_ARG_ta)
        dd <- core_ARG_ta %>%
            group_by(Group) %>%
            summarise(Max = max(Abun))
        test <- data.frame(Group = levels(core_ARG_ta$Group),
                           value.x = ifelse(fit1$p.value < 0.05,c("a","b"),c("a","a")),
                           value.y = dd$Max + max(dd$Max)*0.05)
    }else{
        fit1 <- aov(Abun~Group,data = core_ARG_ta)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,alpah=0.05)
        dd <- core_ARG_ta %>%
            group_by(Group) %>%
            summarise(Max = max(Abun))
        test <- data.frame(Group = levels(core_ARG_ta$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max + max(dd$Max)*0.05)
    }
    
    if (length(levels(core_ARG_ta$Group)) == 2) {
        pvalue <- fit1
    }else{
        pvalue <- tuk1
    }
    
    aa <- max(str_length(levels(group$Group)))
    core_ARG_ta1 <- ggplot(core_ARG_ta,aes(Group,Abun,fill = Group)) + 
        geom_boxplot() +
        labs(y = "Group distance",x = "") +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        scale_fill_manual(values = cbbPalette) + 
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 10),
              axis.text.y=element_text(colour='black',size = 10),
              axis.text.x=element_text(colour = "black",size = 12,face = "bold",
                                       angle = 90,hjust = 1,vjust = 0.5),
              legend.position = "none")
    
    core_ARG_ta2 <- ggplot(core_ARG_ta,aes(Group,Abun,color = Group)) + 
        geom_boxplot(width = 0.5,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        labs(y = "Group distance",x = "") +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black") +
        scale_color_manual(values = cbbPalette) + 
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(color = "black",size = 10),
              axis.text.y=element_text(colour='black',size = 10),
              axis.text.x=element_text(colour = "black",size = 12,
                                       angle = 90,hjust = 1,vjust = 0.5),
              legend.position = "none")
    
    core_ARG_ta3 <- ggplot(core_ARG_ta,aes(Group,Abun,fill = Group)) + 
        geom_violin(scale = "width") +
        labs(y = "Group distance",x = "") +
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        scale_fill_manual(values = cbbPalette) + 
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 10),
              axis.text.y=element_text(colour='black',size = 10),
              axis.text.x=element_text(colour = "black",size = 12,face = "bold",
                                       angle = 90,hjust = 1,vjust = 0.5),
              legend.position = "none")
    result <- list(core_ARG_ta1,core_ARG_ta2,core_ARG_ta3,pvalue)
    return(result)
}