sankey.core.species <- function(ARG_sub,otu3){
    otu3$Phylum <- ifelse(otu3$Phylum == otu3$Class & otu3$Phylum != "Unclassified",
                          paste(otu3$Phylum,"_P",sep = ""),otu3$Phylum)
    otu3$Class <- ifelse(otu3$Class == otu3$Order & otu3$Class != "Unclassified",
                          paste(otu3$Class,"_C",sep = ""),otu3$Class)
    otu3$Order <- ifelse(otu3$Order == otu3$Family & otu3$Order != "Unclassified",
                          paste(otu3$Order,"_O",sep = ""),otu3$Order)
    otu3$Family<- ifelse(otu3$Family == otu3$Genus & otu3$Family != "Unclassified",
                          paste(otu3$Family,"_F",sep = ""),otu3$Family)
    otu3$Genus<- ifelse(otu3$Genus == otu3$Species & otu3$Genus != "Unclassified",
                          paste(otu3$Genus,"_G",sep = ""),otu3$Genus)
    
    ARG_sub$Phylum <- ifelse(ARG_sub$Phylum == ARG_sub$Class & ARG_sub$Phylum != "Unclassified",
                          paste(ARG_sub$Phylum,"_P",sep = ""),ARG_sub$Phylum)
    ARG_sub$Class <- ifelse(ARG_sub$Class == ARG_sub$Order & ARG_sub$Class != "Unclassified",
                         paste(ARG_sub$Class,"_C",sep = ""),ARG_sub$Class)
    ARG_sub$Order <- ifelse(ARG_sub$Order == ARG_sub$Family  & ARG_sub$Order != "Unclassified",
                         paste(ARG_sub$Order,"_O",sep = ""),ARG_sub$Order)
    ARG_sub$Family<- ifelse(ARG_sub$Family == ARG_sub$Genus & ARG_sub$Family != "Unclassified",
                         paste(ARG_sub$Family,"_F",sep = ""),ARG_sub$Family)
    ARG_sub$Genus<- ifelse(ARG_sub$Genus == ARG_sub$Species & ARG_sub$Genus!= "Unclassified",
                        paste(ARG_sub$Genus,"_G",sep = ""),ARG_sub$Genus)
    
    id.phylum <- data.frame(source = rep("Bacteria",length(unique(otu3$Phylum))),
                            target = unique(otu3$Phylum))
    id.class <- data.frame(source = otu3[!duplicated(otu3[,c("Phylum","Class")]),c("Phylum","Class")]$Phylum,
                           target = otu3[!duplicated(otu3[,c("Phylum","Class")]),c("Phylum","Class")]$Class)
    id.order <- data.frame(source = otu3[!duplicated(otu3[,c("Class","Order")]),c("Class","Order")]$Class,
                           target = otu3[!duplicated(otu3[,c("Class","Order")]),c("Class","Order")]$Order)
    id.family <- data.frame(source = otu3[!duplicated(otu3[,c("Order","Family")]),c("Order","Family")]$Order,
                            target = otu3[!duplicated(otu3[,c("Order","Family")]),c("Order","Family")]$Family)
    id.genus <- data.frame(source = otu3[!duplicated(otu3[,c("Family","Genus")]),c("Family","Genus")]$Family,
                           target = otu3[!duplicated(otu3[,c("Family","Genus")]),c("Family","Genus")]$Genus)
    id.species <- data.frame(source = otu3[!duplicated(otu3[,c("Genus","Species")]),c("Genus","Species")]$Genus,
                             target = otu3[!duplicated(otu3[,c("Genus","Species")]),c("Genus","Species")]$Species)
    ARG_sub$ASV_ID <- rownames(ARG_sub)
    ARG_sub <- ARG_sub[,c("ASV_ID",colnames(ARG_sub)[1:(ncol(ARG_sub)-1)])]
    core_ARG_abun <- ARG_sub %>%
        filter(ASV_ID %in% otu3$ID)
    phylum <- aggregate(core_ARG_abun[,2:(ncol(core_ARG_abun)-6)],
                        list(core_ARG_abun$Phylum),sum)
    phylum <- data.frame(target = phylum$Group.1,value = rowMeans(phylum[,2:ncol(phylum)]))
    core_ARG_abun$Class <- ifelse(core_ARG_abun$Phylum == "Unclassified",NA,core_ARG_abun$Class)
    class <- aggregate(core_ARG_abun[,2:(ncol(core_ARG_abun)-6)],
                        list(core_ARG_abun$Class),sum)
    class <- data.frame(target = class$Group.1,value = rowMeans(class[,2:ncol(class)]))
    class <- class[complete.cases(class),]
    core_ARG_abun$Order <- ifelse(core_ARG_abun$Class == "Unclassified",NA,core_ARG_abun$Order)
    core_ARG_abun$Order <- ifelse(is.na(core_ARG_abun$Class),NA,core_ARG_abun$Order)
    order <- aggregate(core_ARG_abun[,2:(ncol(core_ARG_abun)-6)],
                        list(core_ARG_abun$Order),sum)
    order <- data.frame(target = order$Group.1,value = rowMeans(order[,2:ncol(order)]))
    order <- order[complete.cases(order),]
    core_ARG_abun$Family <- ifelse(core_ARG_abun$Order == "Unclassified",NA,core_ARG_abun$Family)
    core_ARG_abun$Family <- ifelse(is.na(core_ARG_abun$Order),NA,core_ARG_abun$Family)
    family <- aggregate(core_ARG_abun[,2:(ncol(core_ARG_abun)-6)],
                        list(core_ARG_abun$Family),sum)
    family <- data.frame(target = family$Group.1,value = rowMeans(family[,2:ncol(family)]))
    family <- family[complete.cases(family),]
    core_ARG_abun$Genus <- ifelse(core_ARG_abun$Family == "Unclassified",NA,core_ARG_abun$Genus)
    core_ARG_abun$Genus <- ifelse(is.na(core_ARG_abun$Family),NA,core_ARG_abun$Genus)
    genus <- aggregate(core_ARG_abun[,2:(ncol(core_ARG_abun)-6)],
                        list(core_ARG_abun$Genus),sum)
    genus <- data.frame(target = genus$Group.1,value = rowMeans(genus[,2:ncol(genus)]))
    genus <- genus[complete.cases(genus),]
    core_ARG_abun$Species <- ifelse(core_ARG_abun$Genus == "Unclassified",NA,core_ARG_abun$Species)
    core_ARG_abun$Species <- ifelse(is.na(core_ARG_abun$Genus),NA,core_ARG_abun$Species)
    species <- aggregate(core_ARG_abun[,2:(ncol(core_ARG_abun)-6)],
                        list(core_ARG_abun$Species),sum)
    species <- data.frame(target = species$Group.1,value = rowMeans(species[,2:ncol(species)]))
    species <- genus[complete.cases(species),]
    phylum <- merge(id.phylum,phylum)
    phylum <- phylum[phylum$source != "Unclassified",]
    phylum <- phylum[phylum$target != "Unclassified",]
    class <- merge(id.class,class)
    class <- class[class$source != "Unclassified",]
    class <- class[class$target != "Unclassified",]
    order <- merge(id.order,order)
    order <- order[order$source != "Unclassified",]
    order <- order[order$target != "Unclassified",]
    family <- merge(id.family,family)
    family <- family[family$source != "Unclassified",]
    family <- family[family$target != "Unclassified",]
    genus <- merge(id.genus,genus)
    genus <- genus[genus$source != "Unclassified",]
    genus <- genus[genus$target != "Unclassified",]
    species <- merge(id.species,species)
    species <- species[species$source != "Unclassified",]
    species <- species[species$target != "Unclassified",]
    plotdata <- data.frame(rbind(phylum,class,order,family,genus,species))
    return(plotdata)
}