pcoa.community2 <- function(ARG_sub,group,species.type){
    colnames(group) <- c("sample","Group")
    pcoa<- pcoa(ARG_sub, correction = "none", rn = NULL)
    PC1 = pcoa$vectors[,1]
    PC2 = pcoa$vectors[,2]
    plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2)
    colnames(plotdata) <-c("sample","PC1","PC2")
    plotdata <- merge(plotdata,group)
    pc1 <-floor(pcoa$values$Relative_eig[1]*100)
    pc2 <-floor(pcoa$values$Relative_eig[2]*100)
    
    ARG_pcoa1<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title=paste("PCoA - ",species.type," communities",sep = "")) + 
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
              legend.title=element_text(size = 14,face = "bold"),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))
    
    ARG_pcoa2<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        stat_ellipse(aes(fill = Group),geom = "polygon",level = 0.95,alpha = 0.3)+
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title=paste("PCoA - ",species.type," communities",sep = "")) + 
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
              legend.title=element_text(size = 14),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5))
    
    ARG_pcoa3<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        geom_label_repel(aes(PC1,PC2,label = sample),fill = "white",color = "black",
                         box.padding = unit(0.3,"lines"),segment.colour = "grey50",
                         label.padding = unit(0.2,"lines"),size = 2) +
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title=paste("PCoA - ",species.type," communities",sep = "")) + 
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
              legend.title=element_text(size = 14,face = "bold"),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))
    
    yf <- plotdata
    yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
    yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
    yd1$Max <- yd1$Max + max(yd1$Max)*0.1
    yd2$Max <- yd2$Max + max(yd2$Max)*0.1
    
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(plotdata$Group)) == 2) {
        fit1 <- t.test(PC1~Group,data = plotdata)
        fit2 <- t.test(PC2~Group,data = plotdata)
        test <- data.frame(Group = yd1$Group,
                           PC1 = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           PC2 = if(fit2$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           yd1 = yd1$Max,yd2 = yd2$Max)
    }else{
        fit1 <- aov(PC1~Group,data = plotdata)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,level = 0.05)
        fit2 <- aov(PC2~Group,data = plotdata)
        tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
        res2 <- cld(tuk2,level = 0.05)
        
        test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                           yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group)
    }
    test$Group <- factor(test$Group,levels = levels(group$Group))
    p1 <- ggplot(plotdata,aes(Group,PC1)) +
        geom_boxplot(aes(fill = Group)) +
        geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
                  size = 7,color = "black",fontface = "bold") +
        coord_flip() +
        scale_fill_manual(values=cbbPalette) +
        theme_bw()+
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_text(colour='black',size=20,face = "bold"),
              axis.text.x=element_blank(),
              legend.position = "none")
    
    p3 <- ggplot(plotdata,aes(Group,PC2)) +
        geom_boxplot(aes(fill = Group)) +
        geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
                  size = 7,color = "black",fontface = "bold") +
        scale_fill_manual(values=cbbPalette) +
        theme_bw()+
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x=element_text(colour='black',size=20,angle = 45,
                                       vjust = 1,hjust = 1,face = "bold"),
              axis.text.y=element_blank(),
              legend.position = "none")
    
    p2<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=8,pch = 21)+
        scale_fill_manual(values=cbbPalette,name = "Group")+
        xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
        ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
        xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
        ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
        theme(text=element_text(size=30))+
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        theme(panel.background = element_rect(fill='white', colour='black'),
              panel.grid=element_blank(), 
              axis.title = element_text(color='black',size=34),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_text(colour='black', size=34,vjust = 7),
              axis.title.y=element_text(colour='black', size=34,vjust = -2),
              axis.text=element_text(colour='black',size=28),
              legend.title=element_text(size = 24,face = "bold"),
              legend.text=element_text(size=20),
              legend.key=element_blank(),legend.position = c(0.88,0.13),
              legend.background = element_rect(colour = "black"),
              legend.key.height=unit(1,"cm")) +
        guides(fill = guide_legend(ncol = 1))
    
    otu.adonis=adonis(ARG_sub~Group,data = plotdata)
    
    p4 <- ggplot(plotdata, aes(PC1, PC2)) +
        geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",otu.adonis$aov.tab$Df[1],
                                                     "\nR2 = ",round(otu.adonis$aov.tab$R2[1],4),
                                                     "\np-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),
                  size = 7) +
        theme_bw() +
        xlab("") + ylab("") +
        theme(panel.grid=element_blank(), 
              axis.title = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank())
    
    p5 <- p1 + p4 + p2 + p3 + 
        plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
    
    result <- list(plotdata,ARG_pcoa1,ARG_pcoa2,ARG_pcoa3,pcoa$vectors,pcoa,p5)
    return(result)
}