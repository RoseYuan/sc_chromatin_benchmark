# source("scripts/feature_engineering/func_signac.R")
source("scripts/feature_engineering/func_archr.R")
source("scripts/feature_engineering/func_snapatac1.R")
source("scripts/feature_engineering/func_aggregation.R")

saveRdsObject <- function(sobj, path) {
	saveRDS(sobj, file = path)
}

saveFeatureMatrix <- function(mobj, path) {
	write.table(mobj, file = path, sep = "\t", quote = FALSE, col.names = FALSE)
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
	ndim <- 50 # use ndim=50 for clustering
	components <- DepthCorComponentsSignac(sobj, corCutOff = 0.75,
									assay_type="all_cell_peaks", n=ndim,
									reduction="lsi_all_cell_peaks")
	# Do clustering
	sobj <- FindNeighbors(object = sobj,
							 reduction = "lsi_all_cell_peaks",
							 dims = components,
							 graph.name = c(paste0("nn_ndim", ndim), paste0("snn_ndim", ndim)))

	library(reticulate)
	use_python("/home/siluo/Software/mambaforge/bin/python")
	# use_python("/home/siluo/softwares/mambaforge-pypy3/envs/sc-chrom-R4/bin/python") # temporary

	sobj <- FindClusters(object = sobj,
							verbose = FALSE,
							algorithm = 4,
							graph.name = paste0("snn_ndim", ndim))
	sobj[["first_round_clusters"]] <- sobj$seurat_clusters

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

runArchR_tiles <- function(fragfiles, output, genome, scale, resolutions, ndim, tileSize=500){
	suppressPackageStartupMessages({
		require(ArchR)
		require(rhdf5)
		require(parallel)
		require(dplyr)
	})

	n_iteration <- length(resolutions)

	addArchRGenome(genome)
	addArchRVerbose(verbose = FALSE)

	fragfile_list <- strsplit(fragfiles, ",")
	# create arrow file
	ArrowFiles <- createArrowFiles(
	  inputFiles = unlist(fragfile_list) ,
	  sampleNames = paste0("CellinFile",seq(1,length(fragfile_list))),
	  minTSS = 0,
	  minFrags = 0,
	  addTileMat = FALSE,
	  addGeneScoreMat = FALSE,
	  force = FALSE
	)

	# create archR project
	proj <- ArchRProject(
	  ArrowFiles = ArrowFiles,
	  outputDirectory = paste0(output, "/proj"),
	  copyArrows = TRUE #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
	)

	# add Tile matrix
	addTileMatrix(
		input = proj,
	    chromSizes = if (inherits(proj, "ArchRProject")) getChromSizes(proj) else NULL,
	    blacklist = if (inherits(proj, "ArchRProject")) getBlacklist(proj) else NULL,
	    tileSize = tileSize,
	    binarize = TRUE,
	    excludeChr = c("chrM", "chrY"),
	    threads = getArchRThreads(),
	    parallelParam = NULL,
	    force = TRUE,
	    logFile = createLogFile("addTileMatrix")
	)

	n <- nCells(input = proj)
	if (n < 10000){
		nsample <- n
	} else {nsample <- 10000}

	proj <- addIterativeLSI(
		ArchRProj = proj,
	    useMatrix = "TileMatrix",
	    name = paste0("IterativeLSI_ndim",ndim),
	    iterations = n_iteration,
	    clusterParams = list( #See Seurat::FindClusters
		resolution = resolutions, #c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
		sampleCells = nsample,
		n.start = 10
		),
	    varFeatures = 25000,
	    dimsToUse = 1:ndim,
	    scaleDims = scale, # by default
	    force = TRUE
	  )

	return(proj)
	}


runArchR_peaks <- function(fragfiles, output, genome, macs2_path, scale, resolutions, ndim, tileSize=500){
	proj <- runArchR_tiles(fragfiles, output, genome, scale, resolutions, ndim, tileSize)

	# Clustering
	proj <- addClusters(
    input = proj,
    reducedDims = paste0("IterativeLSI_ndim",ndim),
    method = "Seurat",
    name = paste0("Clusters_ndim",ndim),
    resolution = tail(resolutions,n=1),
    algorithm = 4,
    dimsToUse = 1:ndim,
    force = TRUE
    )

	# Adding pseudo-bulk
	proj <- addGroupCoverages(ArchRProj = proj, groupBy = paste0("Clusters_ndim",ndim), force = TRUE)

	# Calling peaks using MACS2
	proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = paste0("Clusters_ndim",ndim),
    peakMethod = "Macs2",
    pathToMacs2 = macs2_path
    )

	# Add insertion matrix
	proj <- addPeakMatrix(proj, binarize = FALSE)

	# Redo LSI
	n_iteration <- length(resolutions)

	n <- nCells(input = proj)
	if (n < 10000){
		nsample <- n
	} else {nsample <- 10000}

	proj <- addIterativeLSI(
	ArchRProj = proj,
	useMatrix = "PeakMatrix",
	name = paste0("IterativeLSI_peaks_ndim",ndim),
	iterations = n_iteration,
	clusterParams = list( #See Seurat::FindClusters
	resolution = resolutions,
	sampleCells = nsample,
	n.start = 10,
	algorithm = 4
	),
	varFeatures = 25000,
	dimsToUse = 1:ndim,
	scaleDims = scale,
	force = TRUE)

	return(proj)
}