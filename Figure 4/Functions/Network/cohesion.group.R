cohesion.group <- function(neg.pos.coh,group){
    ff <- neg.pos.coh %>%
        group_by(Group) %>%
        summarise(mean = mean(value),sd = sd(value))
    
    x <- c("a","b")
    y <- c("a","a")
    if (length(levels(neg.pos.coh$Group)) == 2) {
        fit1 <- t.test(value~Group,data = neg.pos.coh)
        dd <- neg.pos.coh %>%
            group_by(Group) %>%
            summarise(Max = max(value))
        test <- data.frame(Group = levels(neg.pos.coh$Group),
                           value.x = if(fit1$p.value < 0.05){
                               x
                           }else{
                               y
                           },
                           value.y = dd$Max*1.05)
    }else{
        fit1 <- aov(value~Group,data = neg.pos.coh)
        tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
        res1 <- cld(tuk1,level = 0.05)
        dd <- neg.pos.coh %>%
            group_by(Group) %>%
            summarise(Max = max(value))
        test <- data.frame(Group = levels(neg.pos.coh$Group),
                           value.x = res1$mcletters$Letters,
                           value.y = dd$Max*1.05)
    }
    
    if (length(levels(neg.pos.coh$Group)) == 2) {
        aa <- fit1
    }else{
        aa <- summary(tuk1)
    }
    
    
    p <- ggplot(ff, aes(x=Group, y=mean, group=Group, fill=Group)) + 
        geom_bar(stat = "identity",width = 0.7)+
        geom_errorbar(aes(ymin = mean-sd,ymax = mean + sd),
                      width = 0.2) + 
        geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                  size = 5,color = "black",fontface = "bold") +
        scale_fill_manual(values = cbbPalette) +
        labs(y="Negative:Positive cohesion") +
        scale_y_continuous(expand = expansion(mult = c(0,0.1))) + 
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text = element_text(size=12,colour = "black"),
              axis.text.x = element_text(angle=45, hjust=1,face = "bold"),
              axis.title.y = element_text(size=14,colour = "black",face = "bold"),
              axis.title.x = element_blank(),
              legend.position = "none")
    
    result <- list(p,aa)
    return(result)
}