# for a given signac object and a clustering result df file, calculate the value of the evaluation metrics
source("scripts/evaluation/lib_metrics.R")

suppressPackageStartupMessages({
library(optparse)
library(rlang)
library(Signac)
library(Seurat)
})

option_list <- list(
    # parameters for preparing clustering
	make_option(c("-i", "--input"), type="character", default=NA, help="input Signac rds file path."),
	make_option(c("-o", "--output"), type="character", default=NA, help="output result rds file path."),
	make_option(c("-c", "--input_clustering_result"), type="character", default=NA, help="input file path for the clustering results."),
	make_option(c("-m", "--output_metric_file"), type="character", default=NA, help="output file path for evaluation metrics.")
)
opt <- parse_args(OptionParser(option_list=option_list))

sobj <- readRDS(opt$input)
true_labels <- sobj$ground_truth
df_clustering <-  read.table(file = opt$input_clustering_result, sep = "\t")
rownames(df_clustering) <- df_clustering$barcode
clustering <- unlist(lapply(Cells(sobj), function(x){df_clustering[x, "clusterings"]}))

result <- evaluation(sobj, true_labels, clustering)
write.table(result$metrics, file = opt$output_metric_file, sep = "\t", quote = FALSE)
saveRDS(result, file=opt$output)