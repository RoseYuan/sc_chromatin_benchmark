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

