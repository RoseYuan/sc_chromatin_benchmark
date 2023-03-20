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
	peaks <- CallPeaks(obj, macs2.path = macs2_path, group.by = group_by)
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
	# annotation_list <- list("hg19" = EnsDb.Hsapiens.v75,
	# 					 "hg38" = EnsDb.Hsapiens.v86,
	# 					 "mm10" = EnsDb.Mmusculus.v79)
	# annotation <- GetGRangesFromEnsDb(ensdb = annotation_list[[genome]])
	# seqlevelsStyle(annotation) <- "UCSC"

	annotation_list <- list("hg19" = "UCSC_annotation_hg19.rds",
						"hg38" = "UCSC_annotation_hg38.rds",
						"mm10" = "UCSC_annotation_mm10.rds")
	annotation <- readRDS(paste0("database/annotation/", annotation_list[[genome]]))  # Temporary
	
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
