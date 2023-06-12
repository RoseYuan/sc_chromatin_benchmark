# input: name of the dataset, embedding file name, clustering parameters
# outputs: a Signac object with clustering results

source("scripts/clustering/lib_clustering.R", chdir=TRUE)
source("scripts/feature_engineering/func_signac.R", chdir=TRUE)

suppressPackageStartupMessages({
library(optparse)
library(rlang)
library(Signac)
library(Seurat)
})


option_list <- list(
    # parameters for preparing clustering
	make_option(c("-i", "--input"), type="character", default=NA, help="input Signac rds file path."),
	make_option(c("-o", "--output"), type="character", default=NA, help="output Signac rds file path. 
    This object will contain ground truth labels, an embedding matrix, and a SNN graph."),
	make_option(c("-t", "--label_table_file"), type="character", default=NA, help="input file path for ground truth label."),
    make_option(c("-b", "--barcode_col"), type="character", default=NA, help="in label_table_file, column name of cell barcode."),
    make_option(c("-l", "--label_col"), type="character", default=NA, help="in label_table_file, column name of cell label."),
    make_option(c("-e", "--embedding_file"), type="character", default=NA, help="input file path for embedding matrix."),
    make_option(c("-z", "--prepare"), action="store_true", default=FALSE, help="if only preparation steps will be performed or not."),
    # parameters for SNN
    make_option(c("-n", "--ndim"), type="double", default=100, help="number of dimensions for the embedding"),
    # parameters for the merged fuzzy simplical set approach
    make_option(c("-u", "--use_umap"), type="character", default=FALSE, help="If use uwot::similarity_graph to build the similarity graph or not. If not, use SNN graph in Seurat."),
    make_option(c("-k", "--k_umap"), type="double", default=10, help="Number of neighbors used in UMAP graph construction."),

    # parameters for clustering
	make_option(c("-r", "--resolution"), type="double", default=0.2, help="Resolution for clustering"),
    make_option(c("-a", "--algorithm"), type="double", default=4, help="Clustering algorithm"),
    make_option(c("-c", "--clustering_output"), type="character", default=NA, help="output file path for clustering result"),
    make_option(c("-s", "--use_seurat"), action="store_true", default=FALSE, help="if use Seurat::Findclusters() to do clustering or not. If not ,use igraph::cluster::leiden()"),
    make_option(c("-q", "--python"), type="character", default=NA, help="python path for r-reticulate.")
)
# -h should be preserved for --help!!!

opt <- parse_args(OptionParser(option_list=option_list))

# check the required arguments are provided
if (is.na(opt$output)) {
  stop("Output file name must be provided. See script usage (--help)")
}
if (is.na(opt$ndim)) {
  stop("Number of dimensions must be provided. See script usage (--help)")
}

if (opt$prepare) {
    sobj <- readRDS(opt$input)
    sobj <- add_labels(sobj, opt$label_table_file, opt$barcode_col, opt$label_col)
    sobj <- add_embedding(sobj, opt$embedding_file)

    graph_name <- paste0("snn_ndim", opt$ndim)
    if (is.null(sobj@graphs[[graph_name]])) {
        name1 <- paste0("nn_ndim", opt$ndim)
        name2 <- paste0("snn_ndim", opt$ndim)
        sobj <- FindNeighbors(object = sobj, 
                                reduction = "learned_embedding", 
                                graph.name = c(paste0("nn_ndim", opt$ndim), paste0("snn_ndim", opt$ndim))
                            )
        # sobj@graphs[[name1]] <- as.Graph(sobj@graphs[[name1]])
        # sobj@graphs[[name2]] <- as.Graph(sobj@graphs[[name2]])
    }
    ndim0 <- dim(sobj@reductions[["learned_embedding"]])[2]
    sobj <- RunUMAP(sobj, 
                reduction = "learned_embedding",
                dims = 1:ndim0)

    saveRDS(sobj, opt$output)
} else {
    
    if (is.element(opt$use_umap, c("true", "True", "TRUE", "T"))) {
	use_umap <- TRUE
    } else {
        use_umap <- FALSE
        }

    sobj <- readRDS(opt$output)
    if (opt$use_seurat) {
        if (!is.na(opt$python)) {
            if (opt$python != "NA"){
                library(reticulate)
                use_python(opt$python)
            }
        }

    if (!use_umap) {
        sobj <- FindClusters(object = sobj, 
                        verbose = FALSE, 
                        algorithm = opt$algorithm,
                        resolution = opt$resolution,
                        graph.name = paste0("snn_ndim", opt$ndim)
                    )
    } else {
        embed <- Embeddings(Reductions(sobj, "learned_embedding"))
        g <- sobj@graphs[[paste0("snn_ndim", opt$ndim)]]
        sim_graph_adj <- uwot::similarity_graph(embed, n_neighbors = opt$k_umap)
        colnames(sim_graph_adj) <- colnames(g)
        rownames(sim_graph_adj) <- rownames(g)
        sobj@graphs[[paste0("umap_graph_k",opt$k_umap)]] <- as.Graph(sim_graph_adj)

        sobj <- FindClusters(object = sobj, 
                verbose = FALSE, 
                algorithm = opt$algorithm,
                resolution = opt$resolution,
                graph.name = paste0("umap_graph_k",opt$k_umap)
            )
    }
        df_label <- data.frame(sobj$seurat_clusters)
    } else {
        # sce <- as.SingleCellExperiment(snare)
        # graph <- scran::buildSNNGraph(x = sce, use.dimred = "learned_embedding")
        # cluster_leiden <- factor(igraph::cluster_leiden(graph, 
        #                                                 objective_function = "CPM", 
        #                                                 resolution_parameter = opt$resolution, 
        #                                                 n_iterations =3)$membership)
        sobj <- PrepareGraph(sobj, reduction="learned_embedding",
						graph.name.ls=c(paste0("nn_ndim", opt$ndim), paste0("snn_ndim", opt$ndim)), 
						igraph.name=paste0("igraph_snn_ndim", opt$ndim))
        graph <- sobj@misc[[paste0("igraph_snn_ndim", opt$ndim)]]
        cluster_leiden <- factor(igraph::cluster_leiden(graph, 
													objective_function = "modularity",
                                        			resolution_parameter = opt$resolution, 
                                        			n_iterations =10)$membership)

        df_label <- data.frame(cluster_leiden)
    }

    colnames(df_label) <- "clusterings"
    df_label$barcode <- Cells(sobj)
    rownames(df_label) <- NULL
    write.table(df_label, file = opt$clustering_output, sep = "\t", quote = FALSE)
}
