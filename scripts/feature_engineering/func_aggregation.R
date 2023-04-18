source("../feature_engineering/func_signac.R", chdir=TRUE)

# aggregate features
# two-pass mode: useful for large datasets, first do feature clustering on a subset of cells,
# then do feature clustering on meta-cells(over-clustered cells based on the meta-features) 

run_aggregation_method <- function(fragfiles, feature_method, dims=seq(2L, 12L), n_meta_features=1000, n_cells=2000, 
norm_method='tfidf', reduce="pca", ...){

    if (feature_method == "signac_cluster") {
        sobj <- runSignac_ByClusterPeaks(fragfiles, ...)
        feature_matrix <- GetAssayData(sobj[["by_cluster_peaks"]], slot = "counts") # feature-by-cell matrix!
        feature_matrix <- t(feature_matrix)
    } else if (feature_method == "signac_all") {
        sobj <- runSignac_AllCellPeaks(fragfiles, ...)
        feature_matrix <- GetAssayData(sobj[["all_cell_peaks"]], slot = "counts") # feature-by-cell matrix!
        feature_matrix <- t(feature_matrix)
    } else if (feature_method == "archr_tile") {
        sobj <- runArchR_tiles(fragfiles, ...)
        se <- getMatrixFromProject(ArchRProj = sobj, useMatrix = "TileMatrix", verbose = FALSE, binarize = TRUE)
        feature_matrix <- se@assays@data@listData$TileMatrix
        feature_matrix <- t(feature_matrix)
    }else if (feature_method == "archr_peak") {
        sobj <- runArchR_peaks(fragfiles, ...)
        se <- getMatrixFromProject(ArchRProj = sobj, useMatrix = "PeakMatrix", verbose = FALSE, binarize = FALSE)
        feature_matrix <- se@assays@data@listData$TileMatrix
        feature_matrix <- t(feature_matrix)
    }else if (feature_method == "snapatac1") {
        sobj <- runSnapATAC1(fragfiles = fragfiles, ...)
        feature_matrix <- sobj@bmat
    }else if (feature_method == "snapatac2") {
        sobj <- 0 # TODO
    }else{stop("Please specify correct feature method!")}

    if (norm_method == "tfidf") {
        norm_function <- Signac::RunTFIDF
    }
    else{stop("Please specify correct normalization method!")}
    if(is.element(tolower(feature_method), c("signac_all", "signac_cluster"))){
        sce <- as.SingleCellExperiment(sobj)
        res <- aggregate_features(feature_matrix=NULL, dims, n_meta_features, n_cells, norm_function, reduce, sce)
        agg_feature_matrix <- res$Fmat
        sobj[[DefaultAssay(sobj)]][["feature_groups"]]<- res$Fgrp
    } else {
        agg_feature_matrix <- aggregate_features(feature_matrix, dims, n_meta_features, n_cells, norm_function, reduce)
        }
    
    return(list(sobj=sobj, Fmat=agg_feature_matrix))
}

aggregate_features <- function(feature_matrix=NULL, dims, n_meta_features, n_cells, norm_function, reduce, sce=NULL){
    require(SingleCellExperiment)
    # feature_matrix: a cell-by-feature matrix
    # first normalize feature_matrix using norm_function, then run PCA (reduce cell dim), then cluster features, 
    # then log-normalize the meta-features, and (optionaly) lastly do dimensional reduction on meta-feature matrix
    if (is.null(feature_matrix) & is.null(sce)){stop("Please specify the feature matrix or sce object as input!")}
    if (is.null(sce)){
        agg_counts <- scDblFinder:::aggregateFeatures(
        t(feature_matrix),
        dims.use = dims,
        k = n_meta_features,
        num_init = 3,
        use.subset = n_cells,
        norm.fn=norm_function, 
        twoPass=TRUE)
    } else {
        sce_x <- scDblFinder:::aggregateFeatures(
        sce,
        dims.use = dims,
        k = n_meta_features,
        num_init = 3,
        use.subset = n_cells,
        norm.fn=norm_function, 
        twoPass=TRUE)

        agg_counts <- counts(sce_x)
        feature_groups <- metadata(sce_x)$featureGroups
    }

    # create sce object
    sce <- SingleCellExperiment(list(counts=agg_counts))
    # normalize the meta-features
    sce <- scuttle::logNormCounts(sce)

    if (reduce == "original") {
        Fmat <- t(as.matrix(logcounts(sce)))
    }else if (reduce == "pca") {
        pca <- scater::runPCA(t(logcounts(sce)), center=TRUE, scale=TRUE, rank=100)
        Fmat <- as.matrix(pca$x)
    }else if (reduce == 'lsi') {
        tf.idf <- RunTFIDF(object=t(as.matrix(logcounts(sce))), method=1)
        agg_lsi <- RunSVD(t(tf.idf), n = 100,nscale.embeddings = TRUE)
        Fmat <- agg_lsi
    }else{stop("Please specify correct dimensional reduction method!")}

    if (is.null(sce)){
        return(Fmat)
        } else {
            return(list(Fmat=Fmat, Fgrp=feature_groups))
            }
}
