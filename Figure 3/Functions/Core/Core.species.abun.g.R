core.species.abun.g <- function(ARG_sub,core_ARG_id,group){
    ARG_sub <- as.data.frame(ARG_sub)
    ARG_sub$ASV_ID <- rownames(ARG_sub)
    ARG_sub <- ARG_sub[,c("ASV_ID",colnames(ARG_sub)[1:(ncol(ARG_sub)-1)])]
    ARG_sub1 <- melt(ARG_sub)
    colnames(group) <- c("variable","Group")
    colnames(ARG_sub1)[3] <- "Abun"
    core_ARG_abun <- ARG_sub1 %>%
        filter(ASV_ID %in% core_ARG_id$ID)
    core_ARG_ta <- core_ARG_abun %>%
        group_by(variable) %>%
        summarise(Abun = sum(Abun))
    core_ARG_ta <- merge(core_ARG_ta,group)
    core_ARG_ta$Group <- factor(core_ARG_ta$Group)
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(core_ARG_ta$Group)) == 2) {
        fit1 <- t.test(Abun~Group,data = core_ARG_ta)
        dd <- core_ARG_ta %>%
            group_by(Group) %>%
            summarise(Max = max(Abun))
        test <- data.frame(Group = levels(core_ARG_ta$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.05)
    }else{
        fit1 <- aov(Abun~Group,data = core_ARG_ta)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,level = 0.05)
        dd <- core_ARG_ta %>%
            group_by(Group) %>%
            summarise(Max = max(Abun))
        test <- data.frame(Group = levels(core_ARG_ta$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.05)
    }
    
    aa <- max(str_length(levels(group$Group)))
    
    if (length(levels(core_ARG_ta$Group)) == 2) {
        bb <- fit1
    }else{
        bb <- tuk1
    }
        
    core_ARG_ta1 <- ggplot(core_ARG_ta,aes(Group,Abun,fill = Group)) + 
        geom_boxplot() +
        labs(y = "Relative abundance (%) of shared species",x = "") +
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
                                       angle = ifelse(aa > 4,45,0),
                                       hjust = ifelse(aa > 4,1,0.5),
                                       vjust = ifelse(aa > 4,1,0.5)),
              legend.position = "none")
    
    core_ARG_ta2 <- ggplot(core_ARG_ta,aes(Group,Abun,color = Group)) + 
        geom_boxplot(width = 0.5,outlier.color = "transparent") +
        geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
        labs(y = "Relative abundance (%) of shared species",x = "") +
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
                                       angle = ifelse(aa > 4,45,0),
                                       hjust = ifelse(aa > 4,1,0.5),
                                       vjust = ifelse(aa > 4,1,0.5)),
              legend.position = "none")
    
    core_ARG_ta3 <- ggplot(core_ARG_ta,aes(Group,Abun,fill = Group)) + 
        geom_violin(scale = "width") +
        labs(y = "Relative abundance (%) of shared species",x = "") +
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
                                       angle = ifelse(aa > 4,45,0),
                                       hjust = ifelse(aa > 4,1,0.5),
                                       vjust = ifelse(aa > 4,1,0.5)),
              legend.position = "none")
    result <- list(core_ARG_ta1,core_ARG_ta2,core_ARG_ta3,bb)
    return(result)
}