library(igraph)

plot_subgraph <- function(sobj, type1, type2, main, graph_name="snn_ndim15"){
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