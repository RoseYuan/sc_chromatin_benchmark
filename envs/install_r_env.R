x <- c(
# Bioconductor
    "mbkmeans",
    "SingleCellExperiment",
    "GenomeInfoDb",
    "EnsDb.Hsapiens.v75",
    "EnsDb.Hsapiens.v86",
    "EnsDb.Mmusculus.v79",
    "scDblFinder",
    "GenomicRanges",
    "ComplexHeatmap",
    "rhdf5",
    "scran",
    "BiocParallel",
    "bluster",
# CRAN
    "cluster",
    "dplyr",
    "ggplot2",
    "ggpubr",
    "ggrastr",
    "Matrix",
    "pracma",
    "scales",
    "matrixStats",
    "patchwork",
    "RColorBrewer",
    "paletteer",
    "pals",
    "viridis",
    "Seurat",
    "SeuratData",
    "tidyr",
    "MESS",
    "aricode",
    "Signac",
    "interp",
    "s2",
    "ranger",
    "xgboost",
    "circlize",
    "igraph",
    "optparse",
    "reshape2",
    "fpc",
    "clv",
    "rlang",
    "reticulate",
# GitHub
    "GreenleafLab/ArchR",
    "immunogenomics/lisi",
    "r3fang/SnapATAC"
)

# TO INSTALL ALL DEPENDENCIES:
if (!require(BiocManager)) 
    install.packages("BiocManager")
for (. in x) 
    if (!require(., character.only = TRUE)) 
        BiocManager::install(., ask = FALSE, update = TRUE)

# TO CAPTURE SESSIO INFO:
for (. in x) {
    . <- gsub(".*/", "", .)
    suppressPackageStartupMessages(
        library(., character.only = TRUE))
}
si <- capture.output(sessionInfo())
writeLines(si, args[[1]])