bac.process <- function(bac,group,tree){
    colnames(group) <- c("variable","Group")
    group <- group[group$variable %in% colnames(bac),]
    taxonomy <- bac %>% 
        separate(taxonomy,c("Domain","Phylum","Class","Order","Family","Genus","Species"),"; ")
    taxonomy$Phylum <- gsub("p__","",taxonomy$Phylum)
    taxonomy$Class <- gsub("c__","",taxonomy$Class)
    taxonomy$Order <- gsub("o__","",taxonomy$Order)
    taxonomy$Family <- gsub("f__","",taxonomy$Family)
    taxonomy$Genus <- gsub("g__","",taxonomy$Genus)
    taxonomy$Species <- gsub("s__","",taxonomy$Species)
    bac <- bac[,colnames(bac) %in% group$variable]
    bac <- as.data.frame(cbind(bac,taxonomy[,(ncol(taxonomy)-6):ncol(taxonomy)]))
    bac <- bac[bac$Domain != "Unassignable",]
    bac <- bac[bac$Domain != "Unclassified",]
    bac <- subset(bac,select = -Domain)
    bac <- bac[complete.cases(bac$Phylum),]
    bac <- bac[rowSums(bac[,1:(ncol(bac)-6)]) > 1,]
    bac1 <- t(bac[,1:(ncol(bac)-6)])
    bac2 <- bac[,1:(ncol(bac)-6)]
    bac2 <- t(t(bac2)/colSums(bac2)*100)
    bac3 <- as.data.frame(cbind(bac2,bac[,(ncol(bac)-5):ncol(bac)]))
    for (i in (ncol(bac3)-5):ncol(bac3)) {
        bac3[grepl("norank",bac3[,i]),i] <- "Unclassified"
        bac3[grepl("uncultured",bac3[,i]),i] <- "Unclassified"
        bac3[grepl("unclassified",bac3[,i]),i] <- "Unclassified"
        bac3[grepl("unidentified",bac3[,i]),i] <- "Unclassified"
    }
    
    #colnames(bac1) <- paste("\'",colnames(bac1),"\'",sep = "")
    tree <- prune.sample(bac1,tree)
    result <- list(bac1,bac2,bac3,tree)
    return(result)
}