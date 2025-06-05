flower.g <- function(ARG_sub,group){
    ARG_sub <- as.data.frame(ARG_sub)
    ARG_sub$ASV_ID <- rownames(ARG_sub)
    ARG_sub <- ARG_sub[,c("ASV_ID",colnames(ARG_sub)[1:(ncol(ARG_sub)-1)])]
    ARG_sub1 <- melt(ARG_sub)
    colnames(group) <- c("variable","Group")
    Group_numb <- length(unique(group$Group))
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
    core_ARG_num <- length(core_ARG_id)
    
    flower_plot<- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col,circle_col) {
        par( bty = 'n', ann = F, xaxt = 'n', yaxt ='n', mar = c(1,1,1,1))
        plot(c(0,10),c(0,10),type='n')
        n  <- length(sample)
        deg <- 360 / n
        res <- lapply(1:n, function(t){
            draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                         y = 5 + sin((start + deg * (t- 1)) * pi / 180),
                         col = ellipse_col[t],
                         border = ellipse_col[t],
                         a = a, b = b, angle = deg * (t- 1))
            text(x = 5 + 2.5 * cos((start + deg * (t -1)) * pi / 180),
                 y = 5 + 2.5 * sin((start + deg * (t -1)) * pi / 180),
                 otu_num[t],cex = 1.5)
            
            if (deg * (t - 1) < 180 && deg *(t - 1) > 0 ) {
                text(x = 5 + 3.3 * cos((start + deg * (t- 1)) * pi / 180),
                     y = 5 + 3.3 * sin((start + deg * (t- 1)) * pi / 180),
                     sample[t],
                     srt = deg * (t - 1) - start,
                     adj = 1,
                     cex = 1.5
                )
            } else {
                text(x = 5 + 3.3 * cos((start + deg * (t- 1)) * pi / 180),
                     y = 5 + 3.3 * sin((start + deg * (t- 1)) * pi / 180),
                     sample[t],
                     srt = deg * (t - 1) + start,
                     adj = 0,
                     cex = 1.5
                )
            }                   
        })
        draw.circle(x = 5, y = 5, r = r, col =circle_col, border = NA)
        text(x = 5, y = 5, paste('Core:', core_otu),cex = 1.5)
    }
    
    return(flower_plot(sample = sample_id_ARG,otu_num = ARG_num,core_otu = core_ARG_num,
                start = 90,a = 0.3,b = 2,r = 1,circle_col = "white",
                ellipse_col = colorRampPalette(brewer.pal(12,"Set3"))(Group_numb)))
}