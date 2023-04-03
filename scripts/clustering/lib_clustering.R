# function to add ground truth to a Signac object, and check the cell id
add_labels <- function(sobj, label_table_file, barcode_col, label_col){
    suppressPackageStartupMessages({
        require(Signac)
        require(Seurat)
    })
    cell_id <- Cells(sobj)
    df_label <- read.table(label_table_file, header = T, sep="\t", comment.char = "")
    rownames(df_label) <- df_label[, barcode_col]
    cell_labels <- df_label[cell_id, label_col]
    Idents(sobj) <- factor(cell_labels)
    sobj$ground_truth <- factor(cell_labels)
    return(sobj)
}

# function to add an embedding object to a Signac object
add_embedding <- function(sobj, embedding_file){
    suppressPackageStartupMessages({
        require(Signac)
        require(Seurat)
    })
    cell_id <- Cells(sobj)
    
    embed <- read.table(embedding_file, header = F, row.names = 1, sep="\t", comment.char = "")

    
    rownames(embed) <- gsub("CellinFile[0-9]*\\+", "", rownames(embed))
    rownames(embed) <- gsub("CellinFile[0-9]*\\#", "", rownames(embed))

    # cell_id <- intersect(rownames(embed),cell_id)
    cell_id <- cell_id[toupper(cell_id) %in% toupper(rownames(embed))]
    # take only the intersection of cells
    sobj <- subset(x=sobj, cells=cell_id)
    
    g <- rep(seq_along(cell_id), sapply(cell_id, length))
    embed_id <- g[match(toupper(cell_id), toupper(rownames(embed)))]
    embed <- embed[embed_id,]
    if (all(toupper(rownames(embed)) == toupper(Cells(sobj)))){
        rownames(embed) <- Cells(sobj)
    } else(stop("All cells in the embedding being added must match the cells in the object!"))

    colnames(embed) <- NULL

    sobj@reductions[["learned_embedding"]] <- CreateDimReducObject(embeddings = as.matrix(embed), key = "LSI_", assay = DefaultAssay(sobj))
    return(sobj)
}


# # function to do clustering using given parameters
# signac_do_clustering <- function(sobj, r, ndim, key, clustering=4){
#     graph_name <- paste0("snn_ndim", ndim)

#     if (is.null(sobj@graphs[[graph_name]])) {
#         sobj <- FindNeighbors(object = sobj, 
#                               reduction = "learned_embedding", 
#                               dims = 1:ndim,
#                               graph.name = c(paste0("nn_ndim", ndim), paste0("snn_ndim", ndim))
#                             )
#     }

#     sobj <- FindClusters(object = sobj, 
#                          verbose = FALSE, 
#                          algorithm = clustering,
#                          resolution = r,
#                          graph.name = paste0("snn_ndim", n)
#                         )

#     sobj[[paste0("clusters_ndim",ndim,"_r",r)]] <- sobj$seurat_clusters
#     return(sobj)
# }


# # function to calculate UMAP from a specific embeddings, visualize it and colored by labels
# sobj <- Signac::RunUMAP(sobj, 
#                 reduction = embedding_name,
#                 dims = 1:ndim)



