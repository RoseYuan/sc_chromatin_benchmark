import snapatac2 as snap
import pandas as pd
import os.path
import numpy as np

def run_snapatac2(fragfiles, output_file, tile_size, genome, distance='jaccard', ndim=50, black_list=None, nfeatures=1000000):
    fragfile_list = fragfiles.split(",")
    dirname = os.path.dirname(output_file)
    h5ad_files = []
    # Load data
    for i, fragfile in enumerate(fragfile_list):
        name = "CellinFile" + str(i+1)
        tmp_file = dirname + name + ".h5ad"
        data = snap.pp.import_data(
            fragfile,
            genome=snap.genome.hg38,
            sorted_by_barcode=False
        )
        # Add tile matrix, select features
        snap.pp.add_tile_matrix(data, bin_size=tile_size, chunk_size=1000, n_jobs=10)
        h5ad_files.append((name, tmp_file))
        data.write(tmp_file)

    data = snap.AnnDataSet(adatas=h5ad_files, filename=output_file)
    data.obs_names = data.obs['sample'].to_numpy() + "+" + np.array(data.obs_names)

    blacklist_database = {"hg19": "database/blacklist/Blacklist/lists/hg19-blacklist.v2.bed.gz",
                          "hg38": "database/blacklist/Blacklist/lists/hg38-blacklist.v2.bed.gz",
                          "mm10": "database/blacklist/Blacklist/lists/mm10-blacklist.v2.bed.gz"}
    black_list_y = blacklist_database[genome]

    if black_list is None:
        black_list = black_list_y
    print("The number of features used: " + nfeatures + ".")
    snap.pp.select_features(data, min_cells=10, most_variable=nfeatures, blacklist=black_list)

    # dimensional reduction
    ## compute similarity matrix
    snap.tl.spectral(data, n_comps=ndim, features='selected', sample_size=1.0, distance_metric=distance,
                     feature_weights='idf')

    df = pd.DataFrame(data.obsm['X_spectral'])
    df.index = data.obs_names

    data.close()
    return df
