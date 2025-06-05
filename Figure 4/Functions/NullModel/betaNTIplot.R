bNTIRCPlot = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,RCb  = RCb,bNTI = bNTI,group  = "Group"){
    
    ps = inputMicro(otu,tax,map,tree,ps,group  = group)
    ps
    
    psrare <- ps
    map = as.data.frame(sample_data(psrare))
    map$ID = row.names(map)
    sample_data(psrare) = map
    
    
    # Get habitat metadata and add it to the βNTI then merge with the RCbray dataset
    eco.meta1=data.frame(sample_data(psrare)) %>%
        dplyr::select(ID, Group) %>%
        rename(Sample_1 = ID, Group_1 = Group)
    
    eco.meta2=data.frame(sample_data(psrare)) %>%
        dplyr::select(ID, Group) %>%
        rename(Sample_2 = ID, Group_2 = Group)
    
    # bNTI 匹配第一列和第二列的分组信息
    bNTI.df = inner_join(bNTI, eco.meta1) %>%
        inner_join(eco.meta2)
    
    
    # 合并两个数据
    colnames(RCb) = c("Sample_2","Sample_1","RCb")
    turnover.df = inner_join(bNTI.df, RCb)
    head(turnover.df)
    dim(turnover.df)
    
    
    #--------------合并文件保存
    # write.csv(turnover.df,"./Result/bNTI//bNTI_RCbray.csv")
    
    
    
    #-----按照分组统计作图
    
    #------------bNIT作图
    dim(bNTI.df)
    within.bNTI.df = bNTI.df %>%
        filter(Group_1 == Group_2) %>%
        mutate(Group = Group_1)
    
    head(within.bNTI.df )
    
    # map$Group
    # ecosystem.conv = data.frame(Group = c("KO", "OE", "WT"), Group2 = c("Cropland", "Old-field", "Forest"))
    # within.bNTI.df = left_join(within.bNTI.df, ecosystem.conv)
    # within.bNTI.df$Group2 = factor(within.bNTI.df$Group2, levels = c("Cropland", "Old-field", "Forest"))
    # within.bNTI.df$Group = factor(within.bNTI.df$Group, levels=c("KO", "OE", "WT"))
    within.bNTI.df$Group <- factor(within.bNTI.df$Group)
    eco.bNTI.plot = ggplot(within.bNTI.df, aes(x=Group, y=bNTI)) +
        geom_boxplot(outlier.shape=1) +
        geom_hline(yintercept = 2, linetype=2, size=0.5) +
        geom_hline(yintercept = -2, linetype=2, size=0.5) +
        labs(y="betaNTI") +
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size=12,colour = "black"),
              axis.text.x = element_text(angle=45, hjust=1),
              axis.title.x = element_blank(),
              axis.title = element_text(size=14,colour = "black"))
    eco.bNTI.plot
    
    
    
    # 现在按照RCbray进行分开标记系统发育过程
    eco.turnover.df = turnover.df %>%
        filter(Group_1 == Group_2) %>%
        mutate(Group = Group_1)
    
    # eco.turnover.df = left_join(eco.turnover.df, ecosystem.conv)
    # eco.turnover.df$Group2 = factor(eco.turnover.df$Group2, levels = c("Cropland", "Old-field", "Forest"))
    # eco.turnover.df$Group = factor(eco.turnover.df$Group, levels=c("KO", "OE", "WT"))
    
    
    head(eco.turnover.df )
    
    
    ## Calculate the relative influence of each process
    eco.turnover.df = eco.turnover.df %>%
        mutate(process = ifelse(abs(bNTI) < 2,
                                ifelse(abs(RCb) < 0.95, "Drift",
                                       ifelse(RCb >= 0.95, "Dispersal limited",
                                              ifelse(RCb <= -0.95, "Homogenizing dispersal", "ERROR"))),
                                ifelse(bNTI >= 2, "Heterogeneous selection",
                                       ifelse(bNTI <= -2, "Homogeneous selection", "ERROR"))))
    
    
    eco.turnover.df$process = factor(eco.turnover.df$process, levels = c("Drift",
                                                                         "Dispersal limited", "Homogenizing dispersal",
                                                                         "Heterogeneous selection", "Homogeneous selection"))
    
    head(eco.turnover.df)
    #------计算每个组的系统发育过程中五个部分分别占有的比例
    sum.eco.turnover.df = eco.turnover.df %>%
        group_by(Group, process) %>%
        dplyr::summarize(n_sites = n(),
                         perc=(n()/45)*100) %>%
        as.data.frame
    
    sum.eco.turnover.df$process
    
    ## Fix the look of some of the variable names for the figure. This is purely asthetic
    # sum.eco.turnover.mod.df = sum.eco.turnover.df %>%
    #   mutate(process = gsub(" S", "\ns", process)) %>%
    #   mutate(process = gsub(" D", "\nd", process)) %>%
    #   mutate(process = gsub("Drift", "Drift alone\n ", process))
    
    # sum.eco.turnover.mod.df$process = factor(sum.eco.turnover.mod.df$process,
    #                                          levels = c("Drift alone\n ",
    #                                                     "Homogenizing\ndispersal",
    #                                                     "Variable\nselection",
    #                                                     "Homogeneous\nselection"))
    
    Pallete <- data.frame(process = c("Drift",
                                      "Dispersal limited", "Homogenizing dispersal",
                                      "Heterogeneous selection", "Homogeneous selection"),
                          values = c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442"))
    Pallete <- Pallete[Pallete$process %in% sum.eco.turnover.df$process,]
    sum.eco.turnover.df$Group <- factor(sum.eco.turnover.df$Group)
    eco.turnover.plot = ggplot(sum.eco.turnover.df, aes(x=Group, y=perc, fill=process)) +
        geom_bar(stat="identity", color="black") +
        scale_fill_manual(values = Pallete$values) +
        labs(y="Percent of site pairs", fill="Process") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text = element_text(size=12,colour = "black"),
              axis.text.x = element_text(angle=45, hjust=1),
              axis.title.y = element_text(size=14,colour = "black"),
              axis.title.x = element_blank(),
              legend.key.size = unit(10, "mm"),
              legend.text = element_text(size=12,colour = "black"),
              legend.title = element_text(size=14,colour = "black"))
    eco.turnover.plot
    
    
    # Merge the plots
    Pallete <- data.frame(process = c("Drift",
                                      "Dispersal limited", "Homogenizing dispersal",
                                      "Heterogeneous selection", "Homogeneous selection"),
                          values = brewer.pal(5,"Set3"))
    Pallete <- Pallete[Pallete$process %in% sum.eco.turnover.df$process,]
    
    eco.plot = cowplot::plot_grid(eco.bNTI.plot, eco.turnover.plot,
                                  rel_widths=c(0.6, 1), labels=c("A", "B"))
    eco.plot
    
    p <- ggpiestats(sum.eco.turnover.df,x=process,y = Group,counts=n_sites,
               results.subtitle = F,
               slice.label = 'percentage',
               perc.k = 0.5,
               label.repel=T,
               legend.title = "Ecological Process",
               proportion.test = F,
               bf.message = TRUE,label.args = list(size = 4,label.r = unit(0.1,"lines"))) + 
        scale_fill_manual(values = Pallete$values) + 
        theme(strip.text = element_text(size = 14,colour = "black"),
              legend.text = element_text(size = 12,colour = "black"),
              legend.title = element_text(size = 14,colour = "black"))
    
    sdratio <- eco.turnover.df %>%
        mutate(process = ifelse(abs(bNTI) < 2,"Stochastic","Deterministic"))
    
    sdratio = sdratio %>%
        group_by(Group, process) %>%
        dplyr::summarize(n_sites = n()) %>%
        as.data.frame
    
    sdratio <- sdratio %>%
        spread(process,n_sites)
    
    sdratio$sdratio <- sdratio$Deterministic/sdratio$Stochastic
    sdratio <- sdratio[,c(1,4)]
    
    p2 <- ggplot(sdratio,aes(Group,sdratio,fill = Group)) + 
        geom_bar(stat = "identity",width = 0.7) + 
        scale_fill_manual(values = cbbPalette) +
        labs(y="Deterministic/Stochastic proces") +
        scale_y_continuous(expand = expansion(mult = c(0,0.1))) + 
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text = element_text(size=12,colour = "black"),
              axis.text.x = element_text(angle=45, hjust=1),
              axis.title.y = element_text(size=14,colour = "black"),
              axis.title.x = element_blank(),
              legend.position = "none")
    
    return(list(eco.bNTI.plot, eco.turnover.plot,eco.plot,p,p2,turnover.df,sdratio))
}