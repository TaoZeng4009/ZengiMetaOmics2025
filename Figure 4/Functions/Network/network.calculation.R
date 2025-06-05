network.calculation <- function(bac2,bac3){
    network_data <- t(bac2)
    colnames(network_data) <- paste(rownames(bac2),"(",bac3$Phylum,")",sep = "")
    network_data <- network_data[,colMeans(network_data) >= 0.1]
    f <- function(x) sum(x==0)
    d <- apply(network_data,2,f)
    network_data1 <- network_data[,d < nrow(network_data)*0.4]
    set.seed(111)
    res <- fast_correlate2(network_data1,method = "spearman",adjust = TRUE)
    
    net <- res$r
    net[net >= 0.6] = 1
    net[net <= -0.6] = 1
    net[net > -0.6 & net < 0.6] = 0
    
    res$r[lower.tri(res$r)] = 0
    res$p[lower.tri(res$p)] = 0
    data.r <- melt(res$r)
    data.p <- melt(res$p)
    result <- cbind(data.r,data.p)
    result <- result[,c(1,2,3,6)]
    colnames(result) <- c("Var1","Var2","R","P")
    
    result.1 <- result[result$R > 0.6,]
    result.2 <- result[result$R < -0.6,]
    result <- rbind(result.1,result.2)
    result <- result[result$P < 0.05,]
    result<- result[result$R < 1,]
    result <- result[result$Var1 != result$Var2,]
    result$R[result$R > 0.6] = 1
    result$R[result$R < -0.6] = -1
    colnames(result) <- c("Source","Target","R","P")
    
    nodes <- res %>%
        as_tbl_graph(abs(r) > 0.6, p < 0.05) %>%
        as_tibble(what = "vertices")
    nodes$Phylum <- gsub(".*\\(","",nodes$name)
    nodes$Phylum <- gsub("\\)","",nodes$Phylum)
    
    net1 <- res %>%
        as_tbl_graph(abs(r) > 0.6, p < 0.05) %>%
        activate("nodes") %>%
        mutate(Degree = centrality_degree(),
               Module = paste("Module",group_infomap()),
               Phylum = nodes$Phylum)
    
    p1 <- ggraph(net1,layout = "stress") +
        geom_edge_fan(color="lightblue",show.legend=FALSE) + 
        geom_node_point(aes(size = Degree,fill=Module),shape=21)+ 
        scale_fill_discrete(name = "Module")+
        scale_edge_width(range=c(0.2,1))+
        theme_graph()
        
    
    p2 <- ggraph(net1,layout = "stress") +
        geom_edge_fan(color="lightblue",show.legend=FALSE) + 
        geom_node_point(aes(size = Degree,fill=Phylum),shape=21)+ 
        scale_fill_discrete(name = "Phylum")+
        scale_edge_width(range=c(0.2,1))+
        theme_graph()
    
    result1 <- list(net,result,network_data1,p1,p2)
    return(result1)
}