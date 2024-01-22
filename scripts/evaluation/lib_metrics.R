suppressPackageStartupMessages({
library(tidyr)
library(Seurat)

library(ggplot2)
library(dplyr)
library(paletteer)
library(ggpubr)
library(reshape2)
library(viridis)
library(patchwork)
library(cluster)
library(fpc)
library(clv)
library(MESS)
})

#######################################
# partition-level metrics
#######################################

compare_clusterings_external <- function(true, pred, metric="ARI") {
    # check if the cell id matches
    if (length(true) != length(pred)) {
        stop("Error! The two partitionings should have the same length!")
    }
    # unnormalized
    if (tolower(metric) == "mi") { # mutual information
        return(infotheo::mutinformation(true, pred))
    }
    if (tolower(metric) == "vi") { # variation of information
        return(clevr::variation_info(true, pred))
    }
    if (tolower(metric) == "purity") { # variation of information
        return(funtimes::purity(true, pred))
    }
    if (tolower(metric) == "vdm|mirkin") { # vdm: Van Dongen Measure; mirkin: Mirkin Metric
        return(mclustcomp::mclustcomp(true, pred, types = tolower(metric))$scores)
    }
    # normalized/adjusted
    if (tolower(metric) == "ari") { # adjusted rand index
        return(aricode::ARI(true, pred))
    }
    if (tolower(metric) == "ami") { # adjusted mutual information
        return(aricode::AMI(true, pred))
    }
    if (tolower(metric) == "homogeneity") {
        return(clevr::homogeneity(true, pred))
    }
}


adjusted_wallance_indices <-function(true=NULL, pred=NULL, contigency_res=NULL){  # cannot be calculated when there's singletons
    
    ## get pairs using C
    ## ensure that values of c1 and c2 are between 0 and n1
    if (is.null(contigency_res)){
      res <- aricode::sortPairs(true, pred)
      n <- length(true)
    } else{
      res <- contigency_res
      n <- sum(res$ni.)
    }
    spairs <- n*(n-1)/2 # N
    stot <- sum(choose(res$nij, 2), na.rm=TRUE) # T
    srow <- sum(choose(res$ni., 2), na.rm=TRUE) # P
    scol <- sum(choose(res$n.j, 2), na.rm=TRUE) # Q
    a <- stot
    b <- srow-stot
    c <- scol-stot
    d <- spairs+stot-srow-scol
    aw <- (a*d-b*c)/((a+b)*(b+d))
    av <- (a*d-b*c)/((a+c)*(c+d))
    ari <- 2*(a*d-b*c)/((a+b)*(b+d)+(a+c)*(c+d))
    
    awi <- list()
    avj <- list()
    for (i in sort(unique(res$pair_c1))){
      idx <- which(res$pair_c1 == i)
      term1 <- spairs * sum(choose(res$nij[idx], 2)) 
      term2 <- choose(sum(res$nij[idx]), 2) * scol
      term3 <- choose(sum(res$nij[idx]), 2) * (spairs - scol)
      awi[i+1] <- (term1 - term2) / term3
    }

    for (j in sort(unique(res$pair_c2))){
      idx <- which(res$pair_c2 == j)
      term1 <- spairs * sum(choose(res$nij[idx], 2)) 
      term2 <- choose(sum(res$nij[idx]), 2) * srow
      term3 <- choose(sum(res$nij[idx]), 2) * (spairs - srow)
      avj[j+1] <- (term1 - term2) / term3
    }
    
    # remove NA introduced by singletons
    aw2 <- mean(unlist(awi[!is.na(awi)]))
    av2 <- mean(unlist(avj[!is.na(avj)]))
    ari2 <- 2*aw2*av2/(aw2+av2)
    return(list("AW"=aw, "AV"=av, "ARI"=ari, "AW2"=aw2, "AV2"=av2, "ARI2"=ari2,"Awi"=awi, "Avj" = avj))
}


#######################################
# embedding-level metrics
#######################################

#' Compute the average k nearest neighbor overlap.
#' 
#' @param nn_id1 A N x k integer matrix of the near neighbour indices.
#' @param nn_idw A N x k integer matrix of the near neighbour indices.
#' @returns A numeric velue of the average knn overlap. 
knn_agreement <- function(nn_id1, nn_id2, k) {
    if (!all(dim(nn_id1) == dim(nn_id2))){
        stop("Please use indice matrix of the same shape!")
    }
    if (!k <= ncol(nn_id1)){
        stop("k out of indices.")
    }
    nn_purity <- vector(mode = "numeric", length = nrow(nn_id1))
    for (i in seq_len(length.out = nrow(x = nn_id1))) {
      nn_purity[i] <- sum(nn_id1[i, ] %in% nn_id2[i, ]) / k
    }
    return(mean(nn_purity))
}

#' For two embedding matrix, compute the average knn overlap, and optionally the area under the knn overlap score using a series of k.
#' for example k_vec <- seq(2,dim(embeddings1)[[1]],100)
knn_agreement_area <- function(embeddings1, embeddings2, k_vec, area=FALSE){
  knn_score <- c()
  k_max <- max(k_vec)
  nn_id1 <- RANN::nn2(data = embeddings1, k = k_max + 1)$nn.idx[, 2:(k_max+1)]
  nn_id2 <- RANN::nn2(data = embeddings2, k = k_max + 1)$nn.idx[, 2:(k_max+1)]
  for (k in k_vec) {
    knn_score <- c(knn_score, knn_agreement(nn_id1[,1:k], nn_id2[,1:k], k))
  }
  if (area) {
    area <- auc(k_vec, knn_score, from = min(k_vec), to = max(k_vec), type = "spline")/(max(k_vec)-min(k_vec)*1)
    return(list("auc"=area, "score"=knn_score))
  }
  return(knn_score)
}

cal_distance <- function(x, metric="Euclidean") {
        if (metric == "Euclidean"){
        dist.matrix <- stats::dist(x)
    }else if(metric == "correlation"){
        dist.matrix <- amap::Dist(x, method = "correlation")
    }else if(metric == "manhattan"){
        dist.matrix <- stats::dist(x,method = "manhattan")
    }else if(metric == "cosine"){
        dist.matrix <- stylo::dist.cosine(x)
        }
}

# for a given distance matrix and given labels, calculate silhouette score, plot the silhouette distribution and average silhouette.
silhouette_result <- function(dist.matrix, labels, title="", my_col=NULL){
    suppressPackageStartupMessages({
        require(cluster)
        require(pals)
        # require(stylo)
    })
    if (is.null(my_col)){
        my_col <- unlist(polychrome())[1:max(as.numeric(as.factor(labels)))]
    }

    sil = silhouette(as.numeric(as.factor(labels)), dist.matrix)
    p <- plot(sil, border=NA, col=my_col, main=title)
    avg_sil <- mean(sil[, 3])
    return(list("fig" = p, "avg" = avg_sil, "sil" = sil))
}


lisi_result <- function(x, meta_data, label_colname, perplexity = 30, nn_eps = 0, print.summary=FALSE, main=''){
    suppressPackageStartupMessages({
        require(lisi)
    })
    res <- compute_lisi(x, meta_data, label_colname, perplexity, nn_eps)
    df_res <- cbind(res[label_colname], meta_data[label_colname])
    colnames(df_res) <- c("cLISI", label_colname)
    df_res <- df_res[order(df_res[,label_colname], -df_res$cLISI), ]
    df_res[,label_colname] <- as.factor(df_res[,label_colname])
    df_res$name <- as.factor(1:nrow(df_res))
    mapping <- aes_string(x="name", y="cLISI", color=label_colname, fill=label_colname)
    options(repr.plot.width=15, repr.plot.height=5)
    p <- ggplot(df_res, mapping) + 
    geom_bar(stat = "identity")+
    labs(y = "cLISI", x = "",
         title = paste0("Clusters cLISI plot, ",
                        main,
                        "\n Average cLISI: ", 
                        round(mean(df_res$cLISI), 2)))
    
    p <- p + theme(axis.text.x = element_blank(), 
                   axis.ticks.x = element_blank())
    # Print summary
    ave <- tapply(df_res$cLISI, df_res[,label_colname], mean)
    n <- tapply(df_res[,label_colname], df_res[,label_colname], length)
    sum <- data.frame(cluster = names(ave), size = n,
                      ave.cLISI = round(ave,2), stringsAsFactors = TRUE)
    if(print.summary) {print(sum)}
    res_avg <- mean(df_res$cLISI)
    return(list(plot=p, avg=res_avg, lisi=df_res$cLISI))
}


#######################################
#graph-level metrics
#######################################

# compute comminity strength/graph connectivity/PWC
# graph: the whole SNN graph; label_idx: the index of the label of the community to compute
pwc <- function(graph, label, label_idx){
  suppressPackageStartupMessages({
    require(igraph)
  })
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
      b <- sum(df_tmp$from %in% v_1 & df_tmp$to %in% v_2) + sum(df_tmp$to %in% v_1 & df_tmp$from %in% v_2)# between cluster edges of node i
      c <- sum(df_tmp$weight[df_tmp$from %in% v_1 & df_tmp$to %in% v_1])
      d <- sum(df_tmp$weight[df_tmp$from %in% v_1 & df_tmp$to %in% v_2]) + sum(df_tmp$weight[df_tmp$to %in% v_1 & df_tmp$from %in% v_2])
      if(b>=a){j1 <- j1+1}
      if(d>=c){w1 <- w1+1}
  }
  return(list(j1=j1, j1_frac=j1/l1, w1=w1, w1_frac=w1/l1))
}


graph_from_sobj <- function(sobj, graph_name, embedding_name="learned_embedding", n_neighbors=20){
  if(graph_name == "umap"){
    snn_g <- graph_from_sobj(sobj, names(sobj)[startsWith(names(sobj), "snn")][1])
    embed <- Embeddings(Reductions(sobj, embedding_name))
    sim_graph_adj <- uwot::similarity_graph(embed, n_neighbors = n_neighbors)
    colnames(sim_graph_adj) <- colnames(snn_g)
    rownames(sim_graph_adj) <- rownames(snn_g)
    umap_g <- igraph::graph_from_adjacency_matrix(adjmatrix = sim_graph_adj, mode = "undirected", weighted = TRUE, add.colnames = TRUE)
    return(umap_g)
  }else{
    gm <- sobj@graphs[[graph_name]]
    attributes(gm)$class <- "dgCMatrix"
    g <- igraph::graph_from_adjacency_matrix(adjmatrix = gm, add.colnames = TRUE, mode = "directed", weighted = TRUE)
    return(g)
  }
}

avg_pwc <- function(sobj, label, embedding_name="learned_embedding", n_neighbors=20, knn=FALSE, umap=FALSE){
  label_ls <- unique(label)
  avg_knn <- NA
  avg_umap <- NA
  df_knn <- NA
  df_umap <- NA

  # compute the community strength for SNN graph
  name <- names(sobj)[startsWith(names(sobj), "snn_ndim")][1]
  snn_g <- graph_from_sobj(sobj, name, embedding_name=embedding_name, n_neighbors=n_neighbors)
  df_snn <- data.frame(cell_type=c(), weak_cells=c(), weak_frac=c())
  for(c in 1:length(label_ls)){
    res <- pwc(snn_g, label, label_idx=c)
    df_snn <- rbind(df_snn, list(cell_type=label_ls[c], weak_cells=res$j1, weak_frac=res$j1_frac))
  }
  # calculate the average
  avg_snn <- mean(df_snn$weak_frac)

  ##################
  # compute the community strength for KNN graph
  if(knn){
    name <- names(sobj)[startsWith(names(sobj), "nn_ndim")][1]
    knn_g <- graph_from_sobj(sobj, name, embedding_name=embedding_name, n_neighbors=n_neighbors)
    df_knn <- data.frame(cell_type=c(), weak_cells=c(), weak_frac=c())
    for(c in 1:length(label_ls)){
      res <- pwc(knn_g, label, label_idx=c)
      df_knn <- rbind(df_knn, list(cell_type=label_ls[c], weak_cells=res$j1, weak_frac=res$j1_frac))
    }
    # calculate the average
    avg_knn <- mean(df_knn$weak_frac)
  }

  ##################
  # compute the community strength for umap similarity graph
  if(umap){
    umap_g <- graph_from_sobj(sobj, "umap", embedding_name=embedding_name, n_neighbors=n_neighbors)
    df_umap <- data.frame(cell_type=c(), weak_cells=c(), weak_frac=c())
    for(c in 1:length(label_ls)){
      res <- pwc(umap_g, label, label_idx=c)
      df_umap <- rbind(df_umap, list(cell_type=label_ls[c], weak_cells=res$j1, weak_frac=res$j1_frac))
    }
    # calculate the average
    avg_umap <- mean(df_umap$weak_frac)
  }

  return(list(avg_knn=avg_knn, avg_snn=avg_snn, avg_umap=avg_umap, df_knn=df_knn, df_snn=df_snn, df_umap=df_umap))
}

#######################################
#Library size effect
#######################################

# graph: adjacency matrix
# score: dataframe, use the column named col
# calculate the geary C index
# adapted from peakVI's script: https://zenodo.org/record/4728534/files/multiomics_latent_spaces.ipynb?download=1
geary_c <- function(score, col, embed=NULL, graph=NULL, k=30){
    if (is.null(embed) & is.null(graph)){stop("Either an input matrix or an input graph is required.")}
    if (is.null(graph)){
        knn_res <- RANN::nn2(embed, k=k)
        dist <- knn_res$nn.dists
        knn_idx <- knn_res$nn.idx
        idx1 <- rep(knn_idx[,1],each=k)
        idx2 <- c(t(knn_idx))
        numer <- (dim(dist)[1] - 1) * sum(dist * ((score[idx1, col] - score[idx2, col]) ** 2))
    } else {
        require(Matrix)
        graph <- as(graph, "TsparseMatrix")
        dist <- graph@x
        idx1 <- graph@Dimnames[[1]][graph@i+1]
        idx2 <- graph@Dimnames[[1]][graph@j+1]
        numer <- (length(dist) - 1) * sum(dist * ((score[idx1, col] - score[idx2, col]) ** 2))
    }
    
    denom <- 2 * sum(dist) * sum((score[,col] - mean(score[,col])) ** 2)
    return(numer / denom)
}

# compute for a latent representation matrix, what is the proportion of latent axis that 
# are significantly partially correlated with counts. Suggestion: use log of counts
significant_latent_frac <- function(embed, counts_vec, ground_truth, p_th=0.05, add_one=FALSE){
    ndim <- dim(embed)[2]
    n <- 0

    # Convert the string vector to numbers based on the mapping
    mapping <- unique(ground_truth)
    numbers <- 1:length(mapping)
    names(numbers) <- mapping
    z <- numbers[ground_truth]

    for(i in 1:ndim){
        res <- ppcor::pcor.test(x=embed[,i], y=log(counts_vec), z=z, method="pearson")
        if (p <= p_th){ n <- n + 1}
    }
    if(add_one){
        n <- n + 1
        ndim <- ndim + 1
        }
    return(list(frac=n/ndim)) 
}

# Calculate the log(Geary C index) between log counts and latent space distance
cal_geary_c <- function(sobj, k=20, embedding_name="learned_embedding", n_neighbors = 20){
  embed <- Embeddings(Reductions(sobj, embedding_name))
  col <- paste0("nCount_", DefaultAssay(object = sobj))
  log_counts <- log(sobj[[col]])
  # calculate Geary C index for KNN graph
  c_knn <- log(geary_c(log_counts, col, embed=embed, graph=NULL, k=k))
  return(list(log_geary_c_knn=c_knn))
}

#######################################
#Evenness
#######################################

normalize_freq <- function(clusterings){
  sample <- as.data.frame(table(clusterings))
  Ntot <- sum(sample$Freq)
  sample$norm_freq <- sample$Freq/Ntot
  return(sample)
}


Hill_diversity <- function(clusterings, q){
  inf_threshold <- 100
  if (q == 0){
    Diversity <- length(unique(clusterings))
  }else{
    sample <- normalize_freq(clusterings)
    pi <- sample$norm_freq  # relative  abundance  of  species
    pi <- pi[pi > 0]
    if (q == 1) {  # exponential of Shanon entropy
      Diversity <- exp(sum(-pi * log(pi)))
    } else if (q >= inf_threshold) {
      Diversity <- 1 / max(pi)
    } else {
      Diversity <- sum(pi^q)^(1 / (1 - q))
    }
  }
  return(Diversity)
}


eveness <- function(clusterings, a=1, b=0){
  Eveness = Hill_diversity(clusterings, a) / Hill_diversity(clusterings, b)
  return(Eveness)
}

#######################################
#The main function to evaluate latent space (cell embedding + graph-based)
#######################################


evaluation_latent <- function(sobj, true_labels, embedding_name="learned_embedding", dist_metric="Euclidean", metrics=NULL){
  embed <- Embeddings(Reductions(sobj, embedding_name))
  if (is.null(metrics)){
      metrics <- c("Silhouette_label", "FNS", "cLISI_label", "log_geary_c_knn", "avg_pwc_snn") 
  }
  df_metric <- data.frame(matrix(ncol = 2, nrow = length(metrics)))
  colnames(df_metric) <- c("metric", "value")
  df_metric$metric <- metrics
  rownames(df_metric) <- df_metric$metric

  # calculating Silhouette score
  dist.matrix <- cal_distance(embed, metric=dist_metric)

  metric <- "Silhouette_label"
  title <- paste0("Silhouette, ", "true labels, ", dist_metric, "distance")
  re2 <- silhouette_result(dist.matrix, true_labels, title=title)
  df_metric[metric, "value"] <- re2$avg
  
  # calculating FNS
  metric <- "FNS"
  df_sil <- data.frame(re2$sil)[, c("cluster", "sil_width")]
  saveRDS(list(sil=df_sil, true_labels=true_labels), file="~/public/SiyuanLuo/projects/rebuttal/number_of_features/df_sil.RDS")
  df_sil$cell_type <- unlist(lapply(df_sil$cluster, function(x){levels(factor(true_labels))[x]}))
  df_sil <- df_sil[order(df_sil[, "cluster"], -df_sil$sil_width), ]
  df_sil$negative_sil <- df_sil$sil_width < 0
  df_sil_neg <- df_sil %>% group_by(cell_type) %>% summarise(negative_count = sum(negative_sil), count=n())
  df_sil_neg$negative_fraction <- df_sil_neg$negative_count/df_sil_neg$count
  df_metric[metric, "value"] <- mean(df_sil_neg$negative_fraction)
    
  # calculating LISI score
  metric <- "cLISI_label"
  title <- paste0("cLISI, ", "true labels")
  re4 <- lisi_result(x=embed, 
                  meta_data=data.frame(clusterings=true_labels), 
                  label_colname="clusterings", 
                  perplexity = 30, 
                  nn_eps = 0, 
                  main=title)
  df_metric[metric, "value"] <- re4$avg

  # calculating Geary C index
  res_geary_c <- cal_geary_c(sobj, k=20, embedding_name=embedding_name, n_neighbors = 20)
  df_metric["log_geary_c_knn", "value"] <- res_geary_c[["log_geary_c_knn"]]

  # calculating the graph connectivity (PWC)
  res_cs <- avg_pwc(sobj, true_labels, n_neighbors=20, embedding_name=embedding_name)
  # df_metric["avg_pwc_knn", "value"] <- res_cs[["avg_knn"]]
  df_metric["avg_pwc_snn", "value"] <- res_cs[["avg_snn"]]
  # df_metric["avg_pwc_umap", "value"] <- res_cs[["avg_umap"]]

  print(df_metric)
  return(list(metrics=df_metric, sil=re2, fns=df_sil_neg, lisi=re4, pwc=res_cs))
}


#######################################
#The main function to evaluate clustering results
#######################################

evaluation_clustering <- function(true_labels, clustering, metrics=NULL){
  if (is.null(metrics)){
      metrics <- c("ARI","AMI","MI","VI")
  }

  df_metric <- data.frame(matrix(ncol = 2, nrow = length(metrics)))
  colnames(df_metric) <- c("metric", "value")
  df_metric$metric <- metrics
  rownames(df_metric) <- df_metric$metric
  for (metric in df_metric$metric) {
    df_metric[metric, "value"] <- compare_clusterings_external(true_labels, clustering, metric = metric)
  }

  # decomposed ARI
  res <- adjusted_wallance_indices(true_labels, clustering)
  df_metric["AW", "value"] <- res$AW
  df_metric["AW", "metric"] <- "AW"

  df_metric["AV", "value"] <- res$AV
  df_metric["AV", "metric"] <- "AV"

  df_metric["AW2", "value"] <- res$AW2
  df_metric["AW2", "metric"] <- "AW2"

  df_metric["AV2", "value"] <- res$AV2
  df_metric["AV2", "metric"] <- "AV2"

  df_metric["ARI2", "value"] <- res$ARI2
  df_metric["ARI2", "metric"] <- "ARI2"

  print(df_metric)
  return(list(metrics=df_metric, awav=res))
}
