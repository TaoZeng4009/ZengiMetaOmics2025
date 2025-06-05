core.species.id.g <- function(ARG_sub,group){
    ARG_sub <- as.data.frame(ARG_sub)
    ARG_sub$ASV_ID <- rownames(ARG_sub)
    ARG_sub <- ARG_sub[,c("ASV_ID",colnames(ARG_sub)[1:(ncol(ARG_sub)-1)])]
    ARG_sub1 <- melt(ARG_sub)
    colnames(group) <- c("variable","Group")
    colnames(ARG_sub1)[3] <- "Abun"
    ARG_sub3 <- merge(ARG_sub1,group)
    ARG_sub4 <- ARG_sub3 %>%
        group_by(ASV_ID,Group) %>%
        summarise(Abun = mean(Abun))
    ARG_sub4 <- spread(ARG_sub4,Group,Abun)
    ARG_sub4 <- ifelse(ARG_sub4[,2:ncol(ARG_sub4)] > 0,ARG_sub4$ASV_ID,NA)
    
    sample_id_ARG <- colnames(ARG_sub4)
    ARG_id<- unique(ARG_sub4[,1])
    a <- c()
    for (i in 1:nrow(ARG_sub4)) {
        b <-length(which(ARG_sub4[i,] != ""))
        a <- c(a,b)
    }
    ARG_id_u <- ARG_sub4[,1]
    ARG_id_u <- ARG_id_u[a == 1]
    ARG_id_u<- ARG_id_u[which(ARG_id_u != '')]
    ARG_id<- ARG_id[ARG_id != '']
    core_ARG_id<- ARG_id
    ARG_num<- length(ARG_id_u)
    
    for(i in 2:ncol(ARG_sub4)) {
        ARG_id <- unique(ARG_sub4[,i])
        ARG_id_u <- ARG_sub4[,i]
        ARG_id_u <- ARG_id_u[a == 1]
        ARG_id_u<- ARG_id_u[which(ARG_id_u != '')]
        ARG_id <- ARG_id[ARG_id != '']
        core_ARG_id <- intersect(core_ARG_id,ARG_id)
        ARG_num <- c(ARG_num, length(ARG_id_u))
    }
    core_ARG_id <- na.omit(core_ARG_id)
    return(core_ARG_id)
}