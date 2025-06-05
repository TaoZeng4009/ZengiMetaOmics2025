diff.alpha.jitter <- function(alpha1,group){
    colnames(group) <- c("ID","Group")
    alpha <- merge(alpha1,group)
    result <- list()
    for (i in 2:ncol(alpha1)) {
        alpha.test <- alpha[,c(i,9)]
        alpha.test <- alpha.test[complete.cases(alpha.test[,1]),]
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
                               value.y = dd$Max*1.03)
        }else{
            fit1 <- aov(Num~Group,data = alpha.test)
            tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
            res1 <- cld(tuk1)
            dd <- alpha.test %>%
                group_by(Group) %>%
                summarise(Max = max(Num))
            test <- data.frame(Group = levels(alpha.test$Group),
                               value.x = res1$mcletters$Letters,
                               value.y = dd$Max*1.03)
        }
        
        dd <- max(str_length(levels(group$Group)))
        
        alpha.test <- ggplot(alpha.test,aes(Group,Num,color = Group)) + 
            geom_boxplot(width = 0.5,outlier.color = "transparent") +
            geom_jitter(width = 0.3,size = 1.5,alpha = 0.5) +
            geom_text(data = test,aes(x = Group,y = value.y,label = value.x),
                      size = 5,color = "black") +
            labs(y = ifelse(i == 2,"Observed species",paste(colnames(alpha1)[i],"index")),
                 x = "") +
            scale_color_manual(values = cbbPalette) +
            theme_bw()+
            theme(panel.grid=element_blank(),
                  axis.ticks.length = unit(0.4,"lines"), 
                  axis.ticks = element_line(color='black'),
                  axis.line = element_line(colour = "black"), 
                  axis.title.x=element_blank(),
                  axis.title.y = element_text(color = "black",size = 16),
                  axis.text.y=element_text(colour='black',size=10),
                  axis.text.x=element_text(colour = "black",size = 10,
                                           angle = ifelse(dd > 4,45,0),
                                           hjust = ifelse(dd > 4,1,0.5),
                                           vjust = ifelse(dd > 4,1,0.5)),
                  legend.position = "none")
        result[[i-1]] <- alpha.test
    }
    return(result)
}