beta.test <- function(data,group){
    
    pcoa<- pcoa(data, correction = "none", rn = NULL)
    PC1 = pcoa$vectors[,1]
    PC2 = pcoa$vectors[,2]
    plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2)
    colnames(plotdata) <-c("variable","PC1","PC2")
    plotdata <- merge(plotdata,group)
    
    ARG.adonis <- adonis2(data~Group,data = plotdata)
    ARG.anosim <- with(plotdata,anosim(data,Group))
    ARG.mrpp <- with(plotdata,mrpp(data,Group))
    result <- list(ARG.adonis,ARG.anosim,ARG.mrpp)
    return(result)
}