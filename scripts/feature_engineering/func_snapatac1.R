source("../feature_engineering/func_signac.R", chdir=TRUE)

snapRbind_new <- function (obj1, obj2) 
{
    if (!is.snap(obj1)) {
        stop(paste("Error @snapRbind: obj1 is not a snap object!", 
            sep = ""))
    }
    if (!is.snap(obj2)) {
        stop(paste("Error @snapRbind: obj2 is not a snap object!", 
            sep = ""))
    }
    barcode1 = paste(obj1@file, obj1@barcode, sep = ".")
    barcode2 = paste(obj2@file, obj2@barcode, sep = ".")
    if (length(unique(c(barcode1, barcode2))) < length(barcode1) + 
        length(barcode2)) {
        stop("Error: @snapRbind: identifcal barcodes found in obj1 and obj2!")
    }
    rm(barcode1, barcode2)
    gc()
    if (nrow(obj1@metaData) > 0 && nrow(obj2@metaData) > 0) {
        metaData = rbind(obj1@metaData, obj2@metaData)
    }
    else {
        metaData = data.frame()
    }
    feature1 = obj1@feature
    feature2 = obj2@feature
    if ((length(feature1) == 0) != (length(feature2) == 0)) {
        stop("different feature found in obj1 and obj2!")
    }
    else {
        if (length(feature1) > 0) {
            if (FALSE %in% (feature1$name == feature2$name)) {
                stop("Error: @snapRbind: different feature found in obj1 and obj2!")
            }
            feature = feature1
        }
        else {
            feature = feature1
        }
    }
    gc()
    peak1 = obj1@peak
    peak2 = obj2@peak
    if ((length(peak1) == 0) != (length(peak2) == 0)) {
        stop("different peak found in obj1 and obj2!")
    }
    else {
        if (length(peak1) > 0) {
            if (FALSE %in% (peak1$name == peak2$name)) {
                stop("Error: @snapRbind: different feature found in obj1 and obj2!")
            }
            peak = peak1
        }
        else {
            peak = peak1
        }
    }
    rm(peak1, peak2)
    gc()
    bmat1 = obj1@bmat
    bmat2 = obj2@bmat
    if ((length(bmat1) == 0) != (length(bmat2) == 0)) {
        stop("bmat has different dimentions in obj1 and obj2!")
    }
    else {
        bmat = rbind(bmat1, bmat2)
    }
    rm(bmat1, bmat2)
    gc()
    gmat1 = obj1@gmat
    gmat2 = obj2@gmat
    if ((length(gmat1) == 0) != (length(gmat2) == 0)) {
        stop("gmat has different dimentions in obj1 and obj2!")
    }
    else {
        gmat = rbind(gmat1, gmat2)
    }
    rm(gmat1, gmat2)
    gc()
    pmat1 = obj1@pmat
    pmat2 = obj2@pmat
    if ((length(pmat1) == 0) != (length(pmat2) == 0)) {
        stop("pmat has different dimentions in obj1 and obj2!")
    }
    else {
        pmat = rbind(pmat1, pmat2)
    }
    rm(pmat1, pmat2)
    gc()
    dmat1 = obj1@smat@dmat
    dmat2 = obj2@smat@dmat
    if ((length(dmat1) == 0) != (length(dmat2) == 0)) {
        stop("dmat has different dimentions in obj1 and obj2!")
    }
    else {
        dmat = rbind(dmat1, dmat2)
    }
    rm(dmat1, dmat2)
    gc()
    res = newSnap()
    res@feature = feature
    res@barcode = c(obj1@barcode, obj2@barcode)
    res@file = c(obj1@file, obj2@file)
    res@sample = c(obj1@sample, obj2@sample)
    res@metaData = metaData
    res@bmat = bmat
    res@pmat = pmat
    res@peak = peak
    res@gmat = gmat
    res@smat@dmat = dmat
    res@smat@sdev = obj1@smat@sdev
    return(res)
}


preprocessingSnapATAC1 <- function(fragfiles, genome, genome_sizefile, binsize, data_dir, py_env=NULL){
    fragfile_list <- unlist(strsplit(fragfiles, ","))
    for (fragfile in fragfile_list) {
        message(paste0("Preprocessing ", fragfile, ".............\n"))
        if(is.null(py_env)){
            bash_command <- paste0('bash scripts/feature_engineering/snapatac1_preprocessing.sh '
            , fragfile, ' ', genome, ' ', genome_sizefile, ' ', binsize, ' ', data_dir)
        }else{
            bash_command <- paste0('conda run --no-capture-output -n ', py_env, 'bash scripts/feature_engineering/snapatac1_preprocessing.sh '
            , fragfile, ' ', genome, ' ', genome_sizefile, ' ', binsize, ' ', data_dir)
        }

        system(bash_command)
    }
}

loadSnapFile <- function(snapfile_list, sample.names){
    suppressPackageStartupMessages({
    require(SnapATAC)})

    snap.files <- unlist(snapfile_list)
    x.sp.ls <- lapply(seq(snap.files), function(i){
        createSnap(
            file=snap.files[i],
            sample=sample.names[i]
        )
    })
    names(x.sp.ls) <- sample.names

    # combine snap objects
    x.sp <- Reduce(snapRbind_new, x.sp.ls)
    x.sp@metaData["sample"] <- x.sp@sample
    return(x.sp)
}

addRobustBinMatrix <- function(x.sp, binsize, black_list, th_outlier=0.001, th_rareness=0.05, nfeatures=NULL){
    # Add cell-by-bin matrix
    x.sp <- addBmatToSnap(x.sp, bin.size=binsize)

    # Matrix binarization
    # remove top th_outlier items in the count matrix and then convert the remaining non-zero values to 1
    x.sp <- makeBinary(x.sp, mat="bmat", outlier.filter=th_outlier)
    print(head(black_list))
    # filter out any bins overlapping with blacklist (bed file) to prevent from potential artifacts
    black_list.gr <- GRanges(black_list[,1], IRanges(black_list[,2], black_list[,3]))
    idy <- queryHits(findOverlaps(x.sp@feature, black_list.gr))
    if (length(idy) > 0) {
        x.sp <- x.sp[,-idy, mat="bmat"]
    }

    # remove unwanted chromasomes
    chr.exclude <- seqlevels(x.sp@feature)[grep("random|chrM|chrU", seqlevels(x.sp@feature))]
    idy <- grep(paste(chr.exclude, collapse="|"), x.sp@feature)
    if (length(idy) > 0) {
        x.sp = x.sp[,-idy, mat="bmat"]
    }

    # The coverage of bins roughly obeys a log normal distribution. Remove top 5% rare bins
    bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
    hist(
        bin.cov[bin.cov > 0], 
        xlab="log10(bin cov)", 
        main="log10(Bin Cov)", 
        col="lightblue", 
        xlim=c(0, 5)
    )
    message(paste0("Number of features for SnapATAC: ", nfeatures, ", using threshold: ", th_rareness, "."))
    
    if(!is.null(nfeatures)){
        th_rareness <- 1 - nfeatures/length(bin.cov[bin.cov > 0])
    }
    
    bin.cutoff = quantile(bin.cov[bin.cov > 0], 1-th_rareness)
    idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
    x.sp = x.sp[, idy, mat="bmat"]
    return(x.sp)
}

runSnapATAC1 <- function(fragfiles, output, genome, scale=TRUE, ndim, genome_sizefile=NULL, black_list=NULL, binsize=5000, py_env=NULL, nfeatures=NULL){
    suppressPackageStartupMessages({
    require(SnapATAC)
    require(GenomicRanges)})
    
    blacklist_database <- list(
        "hg19" = "database/blacklist/Blacklist/lists/hg19-blacklist.v2.bed.gz",
        "hg38" = "database/blacklist/Blacklist/lists/hg38-blacklist.v2.bed.gz",
        "mm10" = "database/blacklist/Blacklist/lists/mm10-blacklist.v2.bed.gz"
        )

	black_list_y <- blacklist_database[[genome]]

    chromsize_database <- list(
        "hg38" = "database/genome_size/hg38.chrom.sizes",
        "hg19" = "database/genome_size/hg19.chrom.sizes",
        "mm10" = "database/genome_size/mm10.chrom.sizes"
    )

    genome_sizefile_y <- chromsize_database[[genome]]

    genome_sizefile <- SetIfNull(genome_sizefile, genome_sizefile_y)
    black_list <- SetIfNull(black_list, black_list_y)
    df_black_list <- read.table(black_list, sep = "\t")

    preprocessingSnapATAC1(fragfiles, genome, genome_sizefile, binsize, output, py_env)

    # load data
    fragfile_list <- unlist(strsplit(fragfiles, ","))
    snapfile_list <- lapply(fragfile_list, function(i){
        # i <- gsub("tsv.gz", "snap", i)
        i <- gsub("gz", "snap", i)
        i <- gsub(dirname(i), output, i)
        # i <- paste0(i, ".snap")
        return(i)
    })

    sample.names <- paste0("CellinFile",seq(1,length(snapfile_list)))
    x.sp <- loadSnapFile(snapfile_list, sample.names)
    
    # Add cell-by-bin matrix and do feature selection
    x.sp <- addRobustBinMatrix(x.sp, binsize, df_black_list, th_outlier=0.001, th_rareness=0.05, nfeatures)

    # run diffusion maps without subsampling
    x.sp = runDiffusionMaps(
        obj= x.sp,
        input.mat="bmat", 
        num.eigs=ndim
    )
    if (scale==TRUE) {
        embed <- x.sp@smat@dmat
        x.sp@smat@dmat <- apply(embed, 2, base::scale)
    }
    return(x.sp)
}

getFeatureMatrixSnapATAC <- function(x.sp, corCutOff= 0.75, n = 100){
    suppressPackageStartupMessages({
    require(Matrix)})
    
    # give the dim.reduct matrix cell names
    embed <- x.sp@smat@dmat
    rownames(embed) <- x.sp@metaData$barcode
    counts <- Matrix::rowSums(x.sp@bmat)

	n <- SetIfNull(x = n, y = ncol(x = embed))
	embed <- embed[, seq_len(length.out = n)]
    components <- DepthCorComponents(embed, counts, corCutOff, n)
    mobj <- embed[, components]
    return(mobj)
}