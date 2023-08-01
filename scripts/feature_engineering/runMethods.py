from func_snapatac2 import *
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fragfile_list", help="Full path to the fragment file. If there're multiple ones, seperating by comma", type=str) #choices=[1,2,3,4]
    parser.add_argument("-o", "--output_h5ad_file", help="Full path to the h5ad file.", type=str)
    parser.add_argument("-f", "--output_mtx_file", help="Full path to the embedding file.", type=str)
    parser.add_argument("-m", "--method", help="Method for feature engineering.", type=str)
    parser.add_argument("-t", "--feature_type", help="Feature type.", type=str, default="default")
    parser.add_argument("-n", "--ndim", help="Number of dimension for the embedding space.", type=int)
    parser.add_argument("-g", "--genome", help="Genome name.", type=str)
    parser.add_argument("-d", "--distance", help="Distance metric to use.", type=str)
    parser.add_argument("-l", "--tile_size", help="Tile/bin size for calculating the count matrix.", type=int)
    args = parser.parse_args()

    fragfiles = args.input_fragfile_list
    output_file = args.output_h5ad_file
    matrix_file = args.output_mtx_file
    tile_size = args.tile_size
    genome = args.genome
    distance = args.distance
    method = args.method
    feature_type = args.feature_type

    ndim = args.ndim


    if method.lower() == "snapatac2":
        mobj = run_snapatac2(fragfiles=fragfiles, 
                             output_file=output_file,
                             tile_size=tile_size,
                             genome=genome,
                             distance=distance,
                             ndim=ndim)

        mobj.to_csv(matrix_file, sep='\t',header=False)