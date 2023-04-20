# getScriptPath <- function() {
# 	    cmd.args <- commandArgs()
#     m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
#         script.dir <- dirname(regmatches(cmd.args, m))
#         if (length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
# 	    if (length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
# 	    return(script.dir)
# }

# script_dir <- getScriptPath()
source("scripts/feature_engineering/feature_engineering.R", chdir=TRUE)

suppressPackageStartupMessages({
library(optparse)
library(rlang)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
})


option_list <- list(
	make_option(c("-m", "--method"), type="character", default=NA, help="feature engineering method to use"),
	make_option(c("-i", "--input"), type="character", default=NA, help="input fragment file path. If there're multiple
	files, use comma to separate them"),
	make_option(c("-o", "--output"), type="character", default=NA, help="output file for RDS object"),
	make_option(c("-f", "--feature_output"), type="character", default=NA, help="output file for feature matrix"),
	# parameters for features
	make_option(c("-d", "--distance"), type="character", default=NA, help="distance metric name"),
	make_option(c("-t", "--feature_type"), type="character", default=NA, help="feature type"),
	make_option(c("-l", "--tile_size"), type="double", default=NA, help="Tile/bin size for ArchR/SnapATAC1,2 method"),
	make_option(c("-r", "--resolutions"), type="character", default="0.2,0.6,0.8", help="Resolution list for iterativeLSI in ArchR method"),
	# parameters for peak calling
	make_option(c("-z", "--macs2_path"), type="character", default=NA, help="path to macs2"), # TODO change the character 
	make_option(c("-a", "--min_width"), type="double", default=0, help="minimal peak width"),
	make_option(c("-b", "--max_width"), type="double", default=Inf, help="maximal peak width"),
	# parameters for preprocessing
	make_option(c("-s", "--scaling"), type="character", default=NA, help="scale the embedding or not"),
	make_option(c("-c", "--cutoff"), type="double", default=0.75, help="cutoff for correlation with depth"),
	make_option(c("-n", "--ndim"), type="double", default=100, help="number of dimensions for the embedding"),
	
	## parameters for aggregation method
	make_option(c("-k", "--n_meta_features"), type="double", default=1000, help="number of meta-features"),
	make_option(c("-u", "--n_cells"), type="double", default=2000, help="number of cells used to do first round feature clustering"),
	make_option(c("-j", "--norm_method"), type="character", default='tfidf', help="name of the normalization method"),
	make_option(c("-e", "--reduce"), type="character", default='pca', help="name of the dimensional reduction method"),
	make_option(c("-q", "--feature_method"), type="character", default=NA, help="name of the feature method"),

	# dataset specific parameters
	make_option(c("-g", "--genome"), type="character", default=NA, help="genome version")
)
# -h should be preserved for --help!!!
# -p reserved for number of processes

methods <- c("signac", "archr", "aggregation", "snapatac1", "snapatac2")
feature_methods <- c("signac_all", "signac_cluster", "archr_tile", "archr_peak", "snapatac1", "snapatac2")
feature_types <- list(
	"signac" = c("all_cell_peaks", "by_cluster_peaks"),
	"archr" = c("peaks", "tiles"),
	"aggregation" = c("default"),
	"snapatac1" = c("default"),
	"snapatac2" = c("default")
	)

opt <- parse_args(OptionParser(option_list=option_list))

# check the required arguments are provided
if (is.na(opt$input)) {
  stop("Input file must be provided. See script usage (--help)")
}
if (is.na(opt$output)) {
  stop("Output file name must be provided. See script usage (--help)")
}
if (is.na(opt$feature_output)) {
  stop("feature output file must be provided. See script usage (--help)")
}
if (is.na(opt$genome)) {
	stop("Genome name must be provided. See script usage (--help)")
}
if (!is.element(tolower(opt$method), methods)) {
	stop("Correct method name must be provided. See script usage (--help)")
}
if (is.element(tolower(opt$method), names(feature_types))) {
	if (!is.element(tolower(opt$feature_type), feature_types[[tolower(opt$method)]])) {
	stop("Correct feature type must be provided. See script usage (--help)")
	}
}
if (is.element(opt$scaling, c("true", "True", "TRUE", "T"))) {
	scaling <- TRUE
} else {
	scaling <- FALSE
	}

check_peak_options <- function() {
	if (is.na(opt$macs2_path)) {
		stop("Path to MACS2 must be provided. See script usage (--help)")
	}
}

check_tilesize_options <- function(){
	if (opt$tile_size==0){
		stop("Please specify tile size!")
	}
}

check_aggregation_options <- function() {
	norm_methods <- c("tfidf")
	reduce_methods <- c("pca", "lsi", "original")

	if (!is.element(tolower(opt$feature_method), feature_methods)){
		stop("Correct feature method name must be provided. See script usage (--help)")
	}
	if (!is.na(opt$norm_method)){
		if (!is.element(opt$norm_method, norm_methods)){
			stop("Correct normalization method must be provided!")
		}
	}
	if (!is.na(opt$reduce)){
		if (!is.element(opt$reduce, reduce_methods)){
			stop("Correct dimensional reduction method must be provided!")
		}
	}
}

if (tolower(opt$method) == "signac") {
	check_peak_options()
	if (opt$feature_type == "all_cell_peaks") {
		sobj <- runSignac_AllCellPeaks(fragfiles=opt$input,
							macs2_path=opt$macs2_path,
							genome=opt$genome,
							scale=scaling,
							min_width=opt$min_width,
							max_width=opt$max_width)
		mobj <- getFeatureMatrixSignac(sobj,
									embedding_name="lsi_all_cell_peaks",
									assay = "all_cell_peaks",
									corCutOff = opt$cutoff,
									n = opt$ndim)
	}else if (opt$feature_type == "by_cluster_peaks") {
	   	sobj <- runSignac_ByClusterPeaks(fragfiles=opt$input,
									 macs2_path=opt$macs2_path,
									genome=opt$genome,
									scale=scaling,
									min_width=opt$min_width,
									max_width=opt$max_width)
		mobj <- getFeatureMatrixSignac(sobj,
								   embedding_name="lsi_by_cluster_peaks",
								   assay = "by_cluster_peaks",
								   corCutOff = opt$cutoff,
								   n = opt$ndim)
	}
}


if (tolower(opt$method) == "archr") {
	resolutions <- as.numeric(unlist(strsplit(opt$resolutions, ",")))
	root <- getwd()
	wd <- dirname(opt$output)
	if (!dir.exists(wd)) {
		dir.create(wd, recursive=TRUE)
		}
	setwd(wd)
	if (opt$feature_type == "tiles") {
		check_tilesize_options()
		sobj <- runArchR_tiles(fragfiles=opt$input,
							   output=dirname(opt$output),  # except for saving to RDS files, save to ArchR project directory
							   genome=opt$genome,
							   scale=scaling,
							   resolutions=resolutions,
							   ndim=opt$ndim,
							   tileSize=opt$tile_size)
		mobj <- getFeatureMatrixArchR(sobj,
									  paste0("IterativeLSI_ndim",opt$ndim),
									  opt$ndim,
									  opt$cutoff)
	}else if (opt$feature_type == "peaks") {
		check_peak_options()
		sobj <- runArchR_peaks(fragfiles=opt$input,
							   output=dirname(opt$output),
							   genome=opt$genome,
                               macs2_path=opt$macs2_path,
							   scale=scaling,
							   resolutions=resolutions,
							   ndim=opt$ndim,
							   tileSize=opt$tile_size)
		mobj <- getFeatureMatrixArchR(sobj,
									  paste0("IterativeLSI_peaks_ndim",opt$ndim),
									  opt$ndim,
									  opt$cutoff)
	}
	saveArchRProject(ArchRProj = sobj, load = FALSE)
	setwd(root)
}

if (tolower(opt$method) == "snapatac1") {
	check_tilesize_options()
	sobj <- runSnapATAC1(fragfiles = opt$input, 
						 output = dirname(opt$output), 
						 genome = opt$genome, 
						 scale = scaling, 
						 ndim = opt$ndim, 
						 binsize=opt$tile_size)
	
	mobj <- getFeatureMatrixSnapATAC(sobj, corCutOff= opt$cutoff, n = opt$ndim) #sed -i 's/+/:/' 100.tsv
}

# When running aggregation method, the parameters for the feature method
#  is specified in the aggregation method section. 
# Putting feature_type empty can
# avoid running a feature method

if (tolower(opt$method) == "aggregation") {  
	check_aggregation_options()
	feature_method <- opt$feature_method
	if (is.element(tolower(feature_method), c("signac_all", "signac_cluster"))) {
		check_peak_options()
		params <- list(macs2_path=opt$macs2_path,
					genome=opt$genome,
					scale=scaling,
					min_width=opt$min_width,
					max_width=opt$max_width)
	}
	if (tolower(feature_method) == "archr_tile") {
		check_tilesize_options()
		resolutions <- as.numeric(unlist(strsplit(opt$resolutions, ",")))
		root <- getwd()
		wd <- dirname(opt$output)
		if (!dir.exists(wd)) {
			dir.create(wd, recursive=TRUE)
			}
		setwd(wd)

		params <- list(output=dirname(opt$output),  # except for saving to RDS files, save to ArchR project directory
					genome=opt$genome,
					scale=scaling,
					resolutions=resolutions,
					ndim=opt$ndim,
					tileSize=opt$tile_size)

		params <- c(list(fragfiles = opt$input, 
					feature_method = tolower(feature_method),
					dims=seq(2L, 12L), 
					n_meta_features=opt$n_meta_features, 
					n_cells=opt$n_cells, 
					norm_method=opt$norm_method, 
					reduce=opt$reduce),
					params)

		result_ls <- do.call(run_aggregation_method, params)

		sobj <- result_ls$sobj
		mobj <- result_ls$Fmat
		saveArchRProject(ArchRProj = sobj, load = FALSE)
		setwd(root)
	}
	if (tolower(feature_method) == "archr_peak") {
		check_peak_options()
		resolutions <- as.numeric(unlist(strsplit(opt$resolutions, ",")))
		root <- getwd()
		wd <- dirname(opt$output)
		if (!dir.exists(wd)) {
			dir.create(wd, recursive=TRUE)
			}
		setwd(wd)

		params <- list(output=dirname(opt$output),
					genome=opt$genome,
					macs2_path=opt$macs2_path,
					scale=scaling,
					resolutions=resolutions,
					ndim=opt$ndim,
					tileSize=opt$tile_size)
		
		params <- c(list(fragfiles = opt$input, 
					feature_method = tolower(feature_method),
					dims=seq(2L, 12L), 
					n_meta_features=opt$n_meta_features, 
					n_cells=opt$n_cells, 
					norm_method=opt$norm_method, 
					reduce=opt$reduce),
					params)

		result_ls <- do.call(run_aggregation_method, params)
		
		sobj <- result_ls$sobj
		mobj <- result_ls$Fmat
		saveArchRProject(ArchRProj = sobj, load = FALSE)
		setwd(root)
	}
	if (tolower(feature_method) == "snapatac1") {
		check_tilesize_options()
		params <- list(output = dirname(opt$output), 
					genome = opt$genome, 
					scale = scaling, 
					ndim = opt$ndim, 
					binsize=opt$tile_size)
	}
	if (tolower(feature_method) == "snapatac2") {
	params <- list()
	}

	if (!is.element(tolower(feature_method), c("archr_tile", "archr_peak"))){
		params <- c(list(fragfiles = opt$input, 
					feature_method = tolower(feature_method),
					dims=seq(2L, 12L), 
					n_meta_features=opt$n_meta_features, 
					n_cells=opt$n_cells, 
					norm_method=opt$norm_method, 
					reduce=opt$reduce),
					params)

		result_ls <- do.call(run_aggregation_method, params)
		sobj <- result_ls$sobj
		mobj <- result_ls$Fmat
	}

}

saveRdsObject(sobj, opt$output)
saveFeatureMatrix(mobj, opt$feature_output)