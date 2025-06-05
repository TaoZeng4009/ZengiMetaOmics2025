
Loading_Color <- function() {
    cbbPalette <- c("#B2182B","#56B4E9","#E69F00","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999",
                "#ADD1E5")
    return(cbbPalette)
}


Loading_Group <- function(file_dir) {
    group <- read.table(file = file_dir,header = TRUE,sep = "\t")
    colnames(group) <- c("variable","Group")
    group$Group <- factor(group$Group,levels = unique(group$Group))
    return(group)
}

Loading_Table <- function(table_dir) {
    bac <- read.table(file = table_dir,header = FALSE,row.names = 1,sep = "\t",quote = "")
    colnames(bac) <- bac[1,]
    bac <- bac[-1,]
    for (i in 1:(ncol(bac)-1)) {
        bac[,i] <- as.numeric(bac[,i])
    }
    return(bac)
}

Loading_Tree <- function(tree_dir) {
    tree <- read.tree(tree_dir)
    return(tree)
}



Data_Process <- function(bac,group,tree){

    Group_numb <- length(unique(group[,2]))
    Sample_numb <- length(unique(group[,1]))
    wid <- ifelse(Group_numb < 4,1,0.8)

    bac <- bac[,c(group$variable,"taxonomy")]
    bac <- bac[bac$taxonomy != "Unculture",]
    group <- group[group$variable %in% colnames(bac),]
    if (grepl("; ",bac$taxonomy[1])) {
        taxonomy <- bac %>% 
            separate(taxonomy,c("Domain","Phylum","Class","Order","Family","Genus","Species"),"; ")
    }else{
        taxonomy <- bac %>% 
            separate(taxonomy,c("Domain","Phylum","Class","Order","Family","Genus","Species"),";")
        
    }
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
        bac3[grepl("metagenome",bac3[,i]),i] <- "Unclassified"
        bac3[is.na(bac3[,i]),i] <- "Unclassified"
        bac3[,i] <- gsub("\\[","",bac3[,i])
        bac3[,i] <- gsub("\\]","",bac3[,i])
    }
    
    #colnames(bac1) <- paste("\'",colnames(bac1),"\'",sep = "")
    tree <- prune.sample(bac1,tree)
    result <- list(bac1,bac2,bac3,tree,Group_numb,Sample_numb,wid)
    return(result)
}



