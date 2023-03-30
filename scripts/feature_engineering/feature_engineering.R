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

