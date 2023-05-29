library(ggplot2)
library(dplyr)
library(paletteer)
library(ggpubr)

library(RColorBrewer)
library(scales)
library(pals)

library(igraph)
library(dplyr)

# method color
my_col_m <- brewer.pal(12, "Paired")
# cell type color
my_col_c <- unlist(polychrome())
my_color <- c("#FFFFFF","#FFFF66","#FFCC33","#F57F17","#EE0000","#990000","#660000")
my_color_2 <- c("#FFFFFF","#0066FF","#0033FF")  # For small values to be obvious from 0
my_color_3 <- c("#FFFFFF","#0099FF","#0066FF","#0033FF","#000099")
my_color_4 <- c("#FFFFFF","#FFCCCC","#FF9999","#FF6666","#FF3333","#CC0033")
my_col_m2 <- brewer.pal(10, "Paired")[c(1,7,3:6,9:10)]
my_col_m3 <- brewer.pal(10, "Paired")[c(2,7,3:6,9:10)]


# compute comminity strength
# graph: the whole SNN graph; label_idx: the index of the label of the community to compute
community_strength <- function(graph, label, label_idx){
    require(igraph)
    V(graph)$label <- label
    label_ls <- unique(label)
    df_graph <- as_data_frame(graph, what = c("edges"))
    
    v_1 <- V(graph)[V(graph)$label == label_ls[label_idx]]
    v_2 <- V(graph)[V(graph)$label != label_ls[label_idx]]

    
    l1 <- length(v_1)
    
    j1 <- 0
    w1 <- 0
    for (i in v_1){
        df_tmp <- df_graph[df_graph$to == i | df_graph$from == i,] # edges of node i
        a <- sum(df_tmp$from %in% v_1 & df_tmp$to %in% v_1) # within cluster edges of node i
        b <- sum(df_tmp$from %in% v_1 & df_tmp$to %in% v_2) # between cluster edges of node i
        c <- sum(df_tmp$weight[df_tmp$from %in% v_1 & df_tmp$to %in% v_1])
        d <- sum(df_tmp$weight[df_tmp$from %in% v_1 & df_tmp$to %in% v_2])
        if(b>=a){j1 <- j1+1}
        if(d>=c){w1 <- w1+1}
    }
    return(list(j1=j1, j1_frac=j1/l1, w1=w1, w1_frac=w1/l1))
}

cross_table_plot <- function(ground_truth, clusterings, a=1.3, b=5.7, c=2, m=0, n=0.2){
    x <- unique(ground_truth)
    y <- as.factor(unique(clusterings))
    data <- expand.grid(X=x, Y=y)
    cross_count <- table(ground_truth, clusterings) 

    # cell count in the cross_count table
    data$Z1 <- apply(data, 1, function(x){cross_count[x[["X"]],as.numeric(x[["Y"]])]})
    # log transform Z1
    data$Z2 <- apply(data, 1, function(x){log(cross_count[x[["X"]],as.numeric(x[["Y"]])]+1)})
    # row normalize of Z1
    data <- data %>% group_by(X) %>% mutate(Z3 = 100*Z1/sum(Z1))

    top_cluster <- data %>% group_by(X) %>% top_n(1, Z1)

    unselected_Y <- setdiff(unique(data$Y), unique(top_cluster$Y))
    top_cluster <- top_cluster[order(top_cluster$X),]
    new_levels <- as.numeric(c(as.character(unique(top_cluster$Y)),unselected_Y))
    data$Y <- factor(data$Y, levels=new_levels)

    res <- adjusted_wallance_indices(ground_truth, clusterings)

    df_awi <- do.call(rbind, Map(data.frame, "awi"=res$Awi))
    df_awi$cell_type <- levels(ground_truth)
    df_count <- data.frame(table(ground_truth))
    df_count$Freq <- round(df_count$Freq / sum(df_count$Freq),3)
    colnames(df_count) <- c("cell_type", "frac")
    df_awi <- merge(df_awi, df_count, by=c("cell_type"))

    df_avj <- do.call(rbind, Map(data.frame, "avj"=res$Avj))
    df_avj$cell_type <- levels(clusterings)
    df_avj$cell_type <- factor(df_avj$cell_type, levels=new_levels)
    df_count <- data.frame(table(clusterings))
    df_count$Freq <- round(df_count$Freq / sum(df_count$Freq),3)
    colnames(df_count) <- c("cell_type", "frac")
    df_avj <- merge(df_avj, df_count, by=c("cell_type"))

    main <- ggplot(data, aes(Y, X, fill= Z3)) + 
    geom_tile(colour="black", size=0.1) + 
    scale_fill_gradientn(colours = my_color, breaks=seq(0,100,10), guide = guide_colourbar()) +
    labs(x="ATAC cluster", y="True class", fill = "Cells per true class %") + 
    theme(legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.key.width= unit(a, 'cm'),
        legend.text=element_text(size=12)) + 
    theme(plot.margin = unit(c(0, 0, 0, 1), "cm"),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 15, margin = margin(5,0,0,0)), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 15, margin = margin(0,0,0,0)), 
        axis.text.y = element_text(size = 12),)

    bp.x <- ggplot(data = df_avj, aes(x = cell_type, y = avj, fill=frac)) + 
    geom_bar(stat = "identity", colour = "black", size=0.1) + 
    scale_fill_gradientn(colours = my_color_4, guide = guide_colourbar(), limits = c(m, n)) +
    theme(
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size = 15, margin = margin(10,5,0,0)),
        legend.position = "none",
        # legend.position = "top",
        panel.background = element_blank()) + 
    labs(y = "Homogeneity", fill="Fraction")+ theme(plot.margin = unit(c(0, 0, 0, b), "cm"))

    bp.y <- ggplot(data = df_awi, aes(x = cell_type, y = awi, fill=frac)) +
    geom_bar(stat = "identity", colour = "black", size=0.1) + coord_flip() + # theme_ipsum() + #theme_gray() +
    scale_fill_gradientn(colours = my_color_4, guide = guide_colourbar(), limits = c(m, n)) +
    theme(axis.title.x = element_text(size = 15, margin = margin(5,5,0,0)), 
            axis.text.x = element_text(size = 12),
            axis.text.y = element_blank(), 
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            # legend.position="none",
            panel.background = element_blank()) + 
    labs(y="Completemess", fill="Fraction")+ theme(plot.margin = unit(c(0, 0, c, 0), "cm"))

    df_hm <- data.frame(cols = numeric(0), value = numeric(0))

    gg_empty <- df_hm %>% 
    ggplot(aes(x = cols, y = value)) +
    geom_blank() +
    theme(axis.text = element_blank(),
            axis.title = element_blank(),
            line = element_blank(),
            panel.background = element_blank()) + 
            geom_text(aes(label = paste0("ARI = ", round(res$ARI,3), "\n", "ARI2 = ", round(res$ARI2,3))), x = 0.5, y = 0.5)

    plot <- ggarrange(
    bp.x, gg_empty, main, bp.y,
    nrow = 2, ncol = 2, widths = c(3, 1), heights = c(1, 3)
    )
    return(plot)
}


plot_subgraph <- function(sobj, type1, type2, main, graph_name="snn_ndim15"){
    require(library(igraph))
    g <- sobj@graphs[[graph_name]]
    attributes(g)$class <- "dgCMatrix"
    g1 <- graph_from_adjacency_matrix(adjmatrix = g, mode = "undirected", weighted = TRUE, add.colnames = TRUE)
    ground_truth <- sobj$ground_truth
    V(g1)$cluster <- ground_truth
    idx <- ground_truth == type1 | ground_truth == type2
    s_v1 <- V(g1)[idx]
    g1_s <- subgraph(g1, s_v1)
    V(g1_s)$color <- V(g1_s)$cluster
    V(g1_s)$color <- gsub(type1,"red", V(g1_s)$color)
    V(g1_s)$color <- gsub(type2,"white", V(g1_s)$color)
    plot(simplify(g1_s), vertex.label=NA, vertex.size=2, main=main) # remove loops and duplicated edges layout=layout.fruchterman.reingold
    legend("topleft", legend=c(type1, type2), pch=21, pt.bg=c("red","white"))
}