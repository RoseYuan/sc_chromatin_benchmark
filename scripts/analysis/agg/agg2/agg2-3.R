suppressPackageStartupMessages({
    devtools::load_all(path="/home/siluo/public/SiyuanLuo/projects/benchmark/scripts/feature_engineering/scFeatAgg")
    require(Signac)
    require(Seurat)
    require(SingleCellExperiment)
    require(mbkmeans)
    library(tidyr)
    library(scran)
    library(BiocParallel)
    library(bluster)
})

aggregate_features <- function(feature_matrix=NULL, ndims, res, norm_function, reduce, sce=NULL){
    require(SingleCellExperiment)
    # require(scDblFinder)
    # feature_matrix: a cell-by-feature matrix
    # first normalize feature_matrix using norm_function, then run PCA (reduce cell dim), then cluster features,
    # then log-normalize the meta-features, and (optionaly) lastly do dimensional reduction on meta-feature matrix
    if (is.null(feature_matrix) & is.null(sce)){stop("Please specify the feature matrix or sce object as input!")}
    if (is.null(sce)){
        agg_counts <- scFeatAgg::aggregateFeatures3(
        t(feature_matrix),
        res = res,
        ndims = ndims,
        #norm.fn = norm_function,
        capEmb = 3)
    } else {
        sce_x <- scFeatAgg::aggregateFeatures3(
        sce,
        res = res,
        ndims = ndims,
        #norm.fn = norm_function,
        capEmb = 3)

        agg_counts <- counts(sce_x)
        feature_groups <- metadata(sce_x)$featureGroups
    }
    gc()
    # create sce object
    sce2 <- SingleCellExperiment(list(counts=agg_counts))
    # normalize the meta-features
    sce2 <- scuttle::logNormCounts(sce2)

    if (reduce == "original") {
        Fmat <- t(as.matrix(logcounts(sce2)))
    }else if (reduce == "pca") {
        pca <- scater::runPCA(t(logcounts(sce2)), center=TRUE, scale=TRUE, rank=100)
        Fmat <- as.matrix(pca$x)
    }else if (reduce == 'lsi') {
        tf.idf <- RunTFIDF(object=t(as.matrix(logcounts(sce2))), method=1)
        agg_lsi <- RunSVD(t(tf.idf), n = 100, scale.embeddings = TRUE)
        Fmat <- agg_lsi
    }else{stop("Please specify correct dimensional reduction method!")}

    if (is.null(sce)){
        return(Fmat)
        } else {
            return(list(Fmat=Fmat, Fgrp=feature_groups))
            }
}

# Buenrostro_2018

sobj_file <- "/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/Buenrostro_2018/Buenrostro_2018/feature_engineering/R/Signac/by_cluster_peaks/0/default/15.RDS"

sobj <- readRDS(sobj_file)
sce <- as.SingleCellExperiment(sobj)

res <- aggregate_features(feature_matrix=NULL, ndims=20, res=5, norm_function=Signac::RunTFIDF, reduce="pca", sce)

embed <- res$Fmat
counts <- Matrix::rowSums(feature_matrix)
embed <- embed[, seq_len(length.out = ndim_feature_method)]
components <- DepthCorComponents(embed, counts, 0.75, ndim_feature_method)
agg_feature_matrix <- embed[, components]

sobj[[DefaultAssay(sobj)]][["feature_groups"]]<- res$Fgrp

saveRDS(sobj, file="agg2-3.RDS")
