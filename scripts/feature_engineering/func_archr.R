getFeatureMatrixArchR <- function(proj, embedding_name, ndim, corCutOff = 0.75) {
	require(ArchR)
	mobj <- getReducedDims(
    ArchRProj = proj,
    reducedDims = embedding_name,
    returnMatrix = TRUE,
    scaleDims = FALSE,  # already scaled/not scaled when the reduced dimension is created
    dimsToUse = 1:ndim,
    corCutOff = corCutOff
  )
	return(mobj)
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
		n.start = 10,
		algorithm = 3
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
    algorithm = 3,
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
	algorithm = 3
	),
	varFeatures = 25000,
	dimsToUse = 1:ndim,
	scaleDims = scale,
	force = TRUE)

	return(proj)
}