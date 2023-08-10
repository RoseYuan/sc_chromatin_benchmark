source("../feature_engineering/func_signac.R", chdir=TRUE)
# aggregate features
# two-pass mode: useful for large datasets, first do feature clustering on a subset of cells,
# then do feature clustering on meta-cells(over-clustered cells based on the meta-features) 

run_aggregation_method <- function(fragfiles, feature_method, dims=seq(2L, 12L), n_meta_features=1000, n_cells=30000, 
norm_method='tfidf', reduce="pca", ndim=100, sobj_file=NULL, ...){

    if (feature_method == "signac_cluster") {
        res <- try(sobj <- readRDS(sobj_file))
        if (class(res) == "try-error"){
            cat("No Signac output provided. Run Signac from the beginning.")
            sobj <- runSignac_ByClusterPeaks(fragfiles, ...)
        }

        feature_matrix <- GetAssayData(sobj[["by_cluster_peaks"]], slot = "counts") # feature-by-cell matrix!
        feature_matrix <- t(feature_matrix)
    } else if (feature_method == "signac_all") {
        sobj <- runSignac_AllCellPeaks(fragfiles, ...)
        feature_matrix <- GetAssayData(sobj[["all_cell_peaks"]], slot = "counts") # the slot "counts" are raw counts, the slot "data" are normalized data (if any normalization has been applied)
        feature_matrix <- t(feature_matrix)
    } else if (feature_method == "archr_tile") {
        sobj <- runArchR_tiles(fragfiles, ndim=ndim, ...)
        se <- getMatrixFromProject(ArchRProj = sobj, useMatrix = "TileMatrix", verbose = FALSE, binarize = TRUE)
        feature_matrix <- se@assays@data@listData$TileMatrix
        feature_matrix <- t(feature_matrix)
    }else if (feature_method == "archr_peak") {
        sobj <- runArchR_peaks(fragfiles, ndim=ndim, ...)
        se <- getMatrixFromProject(ArchRProj = sobj, useMatrix = "PeakMatrix", verbose = FALSE, binarize = FALSE)
        feature_matrix <- se@assays@data@listData$TileMatrix
        feature_matrix <- t(feature_matrix)
    }else if (feature_method == "snapatac1") {
        sobj <- runSnapATAC1(fragfiles = fragfiles, ndim=ndim, ...)
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

        embed <- res$Fmat
        counts <- Matrix::rowSums(feature_matrix)
        embed <- embed[, seq_len(length.out = ndim)]
        components <- DepthCorComponents(embed, counts, 0.75, ndim)
        agg_feature_matrix <- embed[, components]
        
        sobj[[DefaultAssay(sobj)]][["feature_groups"]]<- res$Fgrp
    } else {
        embed <- aggregate_features(feature_matrix=feature_matrix, dims=dims, n_meta_features=n_meta_features, n_cells=n_cells, norm_function=norm_function, reduce=reduce, sce=NULL)
        counts <- Matrix::rowSums(feature_matrix)
        embed <- embed[, seq_len(length.out = ndim)]
        components <- DepthCorComponents(embed, counts, 0.75, ndim)
        agg_feature_matrix <- embed[, components]
        }
    
    return(list(sobj=sobj, Fmat=agg_feature_matrix))
}

aggregate_features <- function(feature_matrix=NULL, dims, n_meta_features, n_cells, norm_function, reduce, sce=NULL){
    suppressPackageStartupMessages({
    require(SingleCellExperiment)
    require(scDblFinder)
    })
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

    # capping
    Fmat_cap <- capEmb(Fmat, nMADs=3, hardCap=TRUE)
    rownames(Fmat_cap) <- rownames(Fmat)

    if (is.null(sce)){
        return(Fmat_cap)
        } else {
            return(list(Fmat=Fmat_cap, Fgrp=feature_groups))
            }
}


capEmb <- function(x, nMADs=3, hardCap=FALSE){
  if(is(x,"SingleCellExperiment")){
    reducedDim(x, "PCA") <- capEmb(reducedDim(x, "PCA"), nMADs=nMADs)
    return(x)
  }
  if(is.null(nMADs) || is.infinite(nMADs)) return(x)
  as.matrix(as.data.frame(lapply(as.data.frame(x), FUN=function(x){
    if(nMADs==0) return(sign(x)*sqrt(abs(x)))
    xm <- median(x)
    xmads <- median(abs(x-xm))
    xmin <- xm-nMADs*xmads
    xmax <- xm+nMADs*xmads
    if(hardCap){
      x[which(x>xmax)] <- xmax
      x[which(x<xmin)] <- xmin
      return(x)
    }
    w <- which(x>xmax)
    x[w] <- xmax+sqrt(x[w]-xmax)
    w <- which(x<xmin)
    x[w] <- xmin-sqrt(xmin-x[w])
    x
  })))
}