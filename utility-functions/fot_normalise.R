# normalize a matrix to a fraction of total in each assay (column)
fot.normalise <- function(x, norm.multiplier = 1000000000, barplot = FALSE){
  data.sum <-   apply(x, 2, function(y){sum(y, na.rm=TRUE)})
  if(barplot){
    barplot(log2(data.sum))  
  }
  
  ##### do ppb normalisation
  x.mat <- as.matrix(x)
  x.mat.ppb <- apply(x.mat, 2, function(i) i/sum(i, na.rm = T) * norm.multiplier )
  x.mat.ppb <- as.data.frame(x.mat.ppb)
  colnames(x.mat.ppb) <- paste("ppb.", colnames(x.mat.ppb), sep = "")
  return(x.mat.ppb)
}
