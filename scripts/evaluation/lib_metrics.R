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


compare_clusterings_external <- function(true, pred, metric="ARI") {
    suppressPackageStartupMessages({
        library(aricode)
        library(clevr)
    })
    # check if the cell id matches
    if (length(true) != length(pred)){
        stop("Error! The two partitionings should have the same length!")
    }else if (toupper(metric)=="ARI"){
        return(ARI(true, pred))
    }else if (toupper(metric)=="AMI"){
        return(AMI(true, pred))
    }else if (toupper(metric)=="HOMOGENEITY"){
        return(homogeneity(true, pred))
    }
}
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

# sil1 <- silhouette_score(x, clusterings, metric="Euclidean", title=paste0("Silhouette_",method,"_ndim",n,"_r",r," ATAC Clusters"))
# sil2 <- silhouette_score(x, truth, metric="Euclidean", title=paste0("Silhouette_",method,"_ndim",n,"_r",r," RNA clusters"))


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

#' @metrics1: metrics that do not need embeddings to calculate
#' @metrics2: metrics that need embeddings to calculate
evaluation <- function(sobj, true_labels, clustering, embedding_name="learned_embedding", dist_metric="Euclidean", metrics1=NULL, metrics2=NULL){
    embed <- Embeddings(Reductions(sobj,embedding_name))
    if (is.null(metrics1)){
        metrics1 <- c("ARI","AMI","Homogeneity")
    }
    if (is.null(metrics2)){
        metrics2 <- c("Silhouette", "Silhouette_label", "cLISI", "cLISI_label")
    }
    df_metric <- data.frame(matrix(ncol = 2, nrow = length(metrics1)+length(metrics2)))
    colnames(df_metric) <- c("metric", "value")
    df_metric$metric <- c(metrics1, metrics2)
    rownames(df_metric) <- df_metric$metric
    for (metric in df_metric$metric) {
        if (metric %in% metrics1) {
            df_metric[metric, "value"] <- compare_clusterings_external(true_labels, clustering, metric = metric)
        }
    }
    # calculating Silhouette score
    dist.matrix <- cal_distance(embed, metric=dist_metric)

    metric <- "Silhouette"
    title <- paste0("Silhouette, ", "ATAC clusters, ", dist_metric, "distance")
    re1 <- silhouette_result(dist.matrix, clustering, title=title)
    df_metric[metric, "value"] <- re1$avg

    metric <- "Silhouette_label"
    title <- paste0("Silhouette, ", "true labels, ", dist_metric, "distance")
    re2 <- silhouette_result(dist.matrix, true_labels, title=title)
    df_metric[metric, "value"] <- re2$avg

    # calculating cLISI score
    metric <- "cLISI"
    title <- paste0("cLISI, ", "ATAC clusters")
    re3 <- lisi_result(x=embed, 
                    meta_data=data.frame(clusterings=clustering), 
                    label_colname="clusterings", 
                    perplexity = 30, 
                    nn_eps = 0,
                    main=title)
    df_metric[metric, "value"] <- re3$avg

    metric <- "cLISI_label"
    title <- paste0("cLISI, ", "true labels")
    re4 <- lisi_result(x=embed, 
                    meta_data=data.frame(clusterings=true_labels), 
                    label_colname="clusterings", 
                    perplexity = 30, 
                    nn_eps = 0, 
                    main=title)
    df_metric[metric, "value"] <- re4$avg
    print(df_metric)
    return(list(metrics=df_metric, sil1=re1, sil2=re2, lisi1=re3, lisi2=re4))
}
