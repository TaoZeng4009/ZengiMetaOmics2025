inputMicro = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group"){
    
    if (is.null(otu)&is.null(tax)&is.null(map)) {
        ps = ps
        map = as.data.frame(sample_data(ps))
        map = map[, group]
        colnames(map) = "Group"
        map$Group = as.factor(map$Group)
        sample_data(ps) = map
        map = NULL
    }
    
    if (is.null(ps) ) {
        
        if (!is.null(otu)) {
            head(otu)
            otu = as.matrix(otu)
            str(otu)
            
            ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE))
            
        }  
        
        if (!is.null(tax) ) {
            head(tax)
            tax = as.matrix(tax)
            # taxa_names(tax)
            x = tax_table(tax)
            ps = merge_phyloseq(ps,x)
            ps
        }
        
        
        if (!is.null(map) ){
            
            map = map[group]
            
            map[,group] = as.factor(map[,group] )
            map$Group 
            z  = sample_data(map)
            ps = merge_phyloseq(ps,z)
            ps
        }
        if (!is.null(tree) ) {
            # #导入进化树
            h = phy_tree(tree)
            ps = merge_phyloseq(ps,h)
            ps
        }
        
        
    }
    return(ps)
    
}