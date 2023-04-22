SetIfNull <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  else{
    return(x)
  }
}

DepthCorComponents <- function(embed, counts, corCutOff, n){
	depth.cor <- as.data.frame(cor(x = embed, y = counts))
	# depth.cor$counts <- depth.cor[, 1]
	# depth.cor$Component <- seq_len(length.out = nrow(x = depth.cor))
	component.keep <- !abs(depth.cor) >= corCutOff
	return(seq_len(length.out = n)[component.keep])
}

DepthCorComponentsSignac <- function(object, corCutOff, assay_type = NULL, n = NULL, reduction = "lsi") {
	require(Signac)
	require(Seurat)
	assay_type <- SetIfNull(x = assay_type, y = DefaultAssay(object = object))
	dr <- object[[reduction]]
	embed <- Embeddings(object = dr)
	counts <- object[[paste0("nCount_", assay_type)]]
	embed <- embed[rownames(x = counts), ]
	n <- SetIfNull(x = n, y = ncol(x = embed))
	embed <- embed[, seq_len(length.out = n)]
	return(DepthCorComponents(embed, counts, corCutOff, n))
}

getFeatureMatrixSignac <- function(sobj, embedding_name, assay = NULL, corCutOff = 0.75, n = 100) {
	require(Signac)
	embed <- Embeddings(Reductions(sobj, embedding_name))

	# assess the correlation between each embedding component and sequencing depth.
	components <- DepthCorComponentsSignac(sobj, corCutOff, assay, n, embedding_name)
	mobj <- embed[, components]
	return(mobj)
}

peakCallingSignac <- function(obj, macs2_path, genome, min_width, max_width, group_by = NULL) {
	#obj can be a fragment objects or Seurat objects
	require(Signac)
	require(Seurat)
	require(GenomeInfoDb)
	require(GenomicRanges)
	
	# get blacklist name and annotation
	blacklist_signac <- list("hg19" = blacklist_hg19,
							 "hg38" = blacklist_hg38_unified,
							 "mm10" = blacklist_mm10,
							 "dm3" = blacklist_dm3,
							 "dm6" = blacklist_dm6,
							 "ce10" = blacklist_ce10,
							 "ce11" = blacklist_ce11)
	blacklist <- blacklist_signac[[genome]]
	# call peaks for all cells using MACS2
	peaks <- CallPeaks(obj, macs2.path = macs2_path, group.by = group_by, verbose = FALSE)
	# remove peaks on nonstandard chromosomes and in genomic blacklist regions
	peaks <- GenomeInfoDb::keepStandardChromosomes(peaks, pruning.mode = "coarse")
	peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)

	# Filter out bad peaks based on length
	peakwidths <- width(peaks)
	peaks <- peaks[peakwidths  < max_width & peakwidths > min_width]
	return(peaks)
}

createSignacObj <- function(frags, peaks, genome, assay_type) {
	require(Signac)
	require(Seurat)
	require(GenomeInfoDb)
	require(EnsDb.Hsapiens.v75)
	require(EnsDb.Hsapiens.v86)
	require(EnsDb.Mmusculus.v79)

	# get gene annotations
	annotation_list <- list("hg19" = EnsDb.Hsapiens.v75,
						 "hg38" = EnsDb.Hsapiens.v86,
						 "mm10" = EnsDb.Mmusculus.v79)

	res <- try(annotation <- GetGRangesFromEnsDb(ensdb = annotation_list[[genome]]))
	if (class(res) != "try-error"){
		res <- try(seqlevelsStyle(annotation) <- "UCSC")
	}
	if (class(res) == "try-error") {
		cat("Caught an error during generating UCSC annotation. Use a previously generated one.\n")
		annotation_list <- list("hg19" = "UCSC_annotation_hg19.rds",
						"hg38" = "UCSC_annotation_hg38.rds",
						"mm10" = "UCSC_annotation_mm10.rds")
		annotation <- readRDS(paste0("database/annotation/", annotation_list[[genome]]))  # Temporary
	}
	
	# Quantify peaks
	counts <- FeatureMatrix(
	  fragments = frags,
	  features = peaks
	)

	# Create Seurat object
	assay <- CreateChromatinAssay(counts,
									  fragments = frags,
									  genome = genome,
									  annotation = annotation)
	sobj <- CreateSeuratObject(assay, assay = assay_type)
	DefaultAssay(sobj) <- assay_type

	return(sobj)
}


runSignac_AllCellPeaks <- function(fragfiles, macs2_path, genome, scale, min_width, max_width) {
	# create fragment objects
	require(Signac)
	require(Seurat)
	# fragfile_list <- strsplit(fragfiles, ",")
	fragfile_list <- unlist(strsplit(fragfiles, ","))
	frag_list <- lapply(fragfile_list, function(i){Signac::CreateFragmentObject(path = i)})
	# call peaks
	peak_list <- lapply(frag_list,
						function(i){peakCallingSignac(i, macs2_path, genome, min_width, max_width, group_by = NULL)})

	combined.peaks <- GenomicRanges::reduce(x = unlist(GenomicRanges::GRangesList(peak_list)))
	sobj_list <- lapply(frag_list,
						function(i){createSignacObj(i, combined.peaks, genome, assay_type = "all_cell_peaks")})

	    if (length(sobj_list) == 1){
		sobj <- sobj_list[[1]]
	}else{
		sobj <- merge(
				  x = sobj_list[[1]],
				  y = sobj_list[2:length(sobj_list)],
				  add.cell.ids = paste0("CellinFile",seq(1,length(sobj_list)))
				)
	}

	sobj <- FindTopFeatures(sobj,
                             min.cutoff = "q5",
                             assay = "all_cell_peaks")
	sobj <- RunTFIDF(sobj,
						  assay = "all_cell_peaks",
						  method = 1)  # computes log(TF×IDF)
	sobj <- RunSVD(sobj,
						n = 100,
						assay = "all_cell_peaks",
						reduction.key = "LSI_",
						reduction.name = "lsi_all_cell_peaks",
						scale.embeddings = scale)
	return(sobj)
}

runSignac_ByClusterPeaks <- function(fragfiles, macs2_path, genome, scale, min_width, max_width){
	sobj <- runSignac_AllCellPeaks(fragfiles, macs2_path, genome, scale, min_width, max_width)
	ndim <- 30 # use ndim=30 for clustering
	components <- DepthCorComponentsSignac(sobj, corCutOff = 0.75,
									assay_type="all_cell_peaks", n=ndim,
									reduction="lsi_all_cell_peaks")
	# Do clustering

	##################################
	# library(reticulate)
	# use_python("/home/siluo/Software/mambaforge/bin/python")
	# use_python("/home/siluo/softwares/mambaforge-pypy3/envs/sc-chrom-R4/bin/python") # temporary
    # use_python("/home/siluo/softwares/mambaforge-pypy3/bin/python")
	sobj <- FindNeighbors(object = sobj,
							 reduction = "lsi_all_cell_peaks",
							 dims = components,
							 graph.name = c(paste0("nn_ndim", ndim), paste0("snn_ndim", ndim)))
	
	sobj <- FindClusters(object = sobj,
							verbose = FALSE,
							algorithm = 3,
							graph.name = paste0("snn_ndim", ndim))

	sobj[["first_round_clusters"]] <- sobj$seurat_clusters

	##################################

	# sce <- as.SingleCellExperiment(sobj)
	# graph <- scran::buildSNNGraph(x = sce, use.dimred = "LSI_ALL_CELL_PEAKS", k=20,
	# 							  type = "jaccard")

	##################################

	# sobj <- PrepareGraph(sobj, reduction="lsi_all_cell_peaks",
	# 				components=components, 
	# 				graph.name.ls=c(paste0("nn_ndim", ndim), paste0("snn_ndim", ndim)), 
	# 				igraph.name=paste0("igraph_snn_ndim", ndim))

	# graph <- sobj@misc[[paste0("igraph_snn_ndim", ndim)]]
    # cluster_leiden <- factor(igraph::cluster_leiden(graph, 
	# 												objective_function = "modularity",
    #                                     			resolution_parameter = 0.8, 
    #                                     			n_iterations =10)$membership)

	# sobj[["first_round_clusters"]] <- cluster_leiden

	##################################

	# peak calling
	peaks <- peakCallingSignac(sobj, macs2_path, genome, min_width, max_width, group_by = "first_round_clusters")
	# get fragment objects
	frags <- Signac::Fragments(sobj)
	sobj <- createSignacObj(frags, peaks, genome, assay_type = "by_cluster_peaks")

	sobj <- FindTopFeatures(sobj,
                             min.cutoff = "q5",
                             assay = "by_cluster_peaks")
	sobj <- RunTFIDF(sobj,
					 method = 1,
					 assay = "by_cluster_peaks")
	sobj <- RunSVD(sobj,
					n = 100,
					assay = "by_cluster_peaks",
					reduction.key = "LSI_",
					reduction.name = "lsi_by_cluster_peaks",
					scale.embeddings = scale)

	return(sobj)
}

#' Create igraph-compatible graph and save in Seurat object
#'
#' adapted from https://github.com/joshpeters/westerlund/blob/master/R/functions.R
PrepareGraph <- function(sobj, reduction, graph.name.ls, igraph.name, components=NULL) {
	components <- SetIfNull(components, 1:ncol(sobj@reductions[[reduction]]))
	stopifnot(ncol(sobj@reductions[[reduction]]) >= max(components))
	sobj <- Seurat::FindNeighbors(object = sobj, 
									dims = components,
									reduction = reduction, 
									graph.name = graph.name.ls)

	g <- sobj@graphs[[graph.name.ls[2]]]
	if(class(g)!="dgCMatrix"){
		attributes(g)[[1]] <- NULL
	}
	attributes(g)$class <- "dgCMatrix"
	#adj_matrix <- Matrix::Matrix(as.matrix(object@graphs[[graph.name]]), sparse = TRUE)
	g <- igraph::graph_from_adjacency_matrix(adjmatrix = g, mode = "undirected", weighted = TRUE, add.colnames = TRUE)
	sobj@misc[[igraph.name]] <- g
	return(sobj)
}