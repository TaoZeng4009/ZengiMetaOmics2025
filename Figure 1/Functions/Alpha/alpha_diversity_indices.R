Alpha_diversity_index <- function(x, tree = NULL, base = exp(1)) {
    est <- estimateR(x)
    Obs <-  est[1, ]
    Shannon <- diversity(x, index = 'shannon', base = base)
    Simpson <- diversity(x, index = 'simpson')
    Pielou <- Shannon / log(Obs, base)
    goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
    
    result <- rbind(est, Shannon, Simpson,
                    Pielou, goods_coverage)
    if (!is.null(tree)) {
        Pd <- pd(x, tree, include.root = FALSE)[1]
        Pd <- t(Pd)
        result <- rbind(result, Pd)
    }
    result <- as.data.frame(t(result))
    return(result)
}