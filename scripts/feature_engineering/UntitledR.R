SetIfNull <- function(x, y){
  if(is.null(x)){
    return(y)
  }
  else{
    return(x)
  }
}

DepthCorComponentsSignac <- function(object, corCutOff, assay = NULL, reduction = 'lsi', n = 100) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  dr <- object[[reduction]]
  embed <- Embeddings(object = dr)
  counts <- object[[paste0("nCount_", assay)]]
  embed <- embed[rownames(x = counts), ]
  n <- SetIfNull(x = n, y = ncol(x = embed))
  embed <- embed[, seq_len(length.out = n)]
  depth.cor <- as.data.frame(cor(x = embed, y = counts))
  depth.cor$counts <- depth.cor[, 1]
  depth.cor$Component <- seq_len(length.out = nrow(x = depth.cor))
  
  component.keep <- !abs(depth.cor) >=corCutOff
  return(seq_len(length.out = n)[component.keep])
}