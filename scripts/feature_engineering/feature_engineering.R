# source("scripts/feature_engineering/func_signac.R")
source("../feature_engineering/func_archr.R", chdir=TRUE)
source("../feature_engineering/func_snapatac1.R", chdir=TRUE)
source("../feature_engineering/func_aggregation.R", chdir=TRUE)

saveRdsObject <- function(sobj, path) {
	saveRDS(sobj, file = path)
}

saveFeatureMatrix <- function(mobj, path) {
	write.table(mobj, file = path, sep = "\t", quote = FALSE, col.names = FALSE)
}

