# candidate1
# Signac_all_cell_peaks, n=15,r=0.25 and Signac_by_cluster_peaks, n=15,r=0.25

# ArchR_peaks, n=15, r=0.15, ArchR_tiles, n=15, r=0.15

# SnapATAC1, n=15, r=0.3
library(reticulate)
library(Signac)
library(Seurat)

lm <- c("Signac_all_cell_peaks", "Signac_by_cluster_peaks", "ArchR_peaks", "ArchR_tiles", "SnapATAC1")
r <- c(0.25,0.25,0.15,0.15,0.3)
python_path <- "/home/siluo/softwares/mambaforge-pypy3/envs/r-reticulate/bin/python"
n <- 15
am <- 4
snn_files <- c("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/candidate1/candidate1/clustering/Signac/all_cell_peaks/0/default/15/sobj_SNN.RDS",
               "/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/candidate1/candidate1/clustering/Signac/by_cluster_peaks/0/default/15/sobj_SNN.RDS",
               "/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/candidate1/candidate1/clustering/ArchR/peaks/500/default/15/sobj_SNN.RDS",
               "/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/candidate1/candidate1/clustering/ArchR/tiles/500/default/15/sobj_SNN.RDS",
               "/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/candidate1/candidate1/clustering/SnapATAC1/default/5000/default/15/sobj_SNN.RDS")


df_meta <- data.frame(long_method=lm,
                      resolution=r, 
                      output=snn_files)

use_python(python_path)

seed_ls <- c(197,123,5,2,42)

for(i in 1:(dim(df_meta)[1]-1)) {
    sobj <- readRDS(df_meta$output[i])
    long_method <- df_meta$long_method[i]
    r <- df_meta$resolution[i]
    for (seed in seed_ls){
        sobj <- FindClusters(object = sobj, 
                    verbose = FALSE, 
                    algorithm = am,
                    resolution = r,
                    graph.name = paste0("snn_ndim", n),
                    random.seed = seed
                    )
        sobj[[paste0(long_method, "_ndim", n, "_r", r ,"_seed", seed)]] <- sobj$seurat_clusters
    }
    saveRDS(sobj, paste0("candidate1/", long_method, "_ndim", n, "_r", r, "seeds.RDS"))
}
