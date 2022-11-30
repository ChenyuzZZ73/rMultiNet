#' @title     A function for tensor decomposition
#' @name      SpecClustering
#' @param     matrx  an n by n matrix
#
#' @param       ncol number of columns of the output matrix U
#' @import      Matrix
#' @returns     U  n by rank numeric matrix that contains the rank tops eigenvectors of Laplacian matrix as column
##' @export


SpecClustering <- function(tnsr, rank, embedding_type = "Layer")
{
  if(embedding_type == "Layer")
  {
    k = k_unfold(tnsr, m=3)
    matrx = k@data
  }
  if(embedding_type == "Nodes")
  {
    matrx = modeSum(tnsr@data, m = 3, drop = TRUE)
  }
  #calculate degree matrix
  diag <- apply(matrx, 1 , sum)
  l <- length(diag)
  D <- sparseMatrix(1:l,1:l,x=diag)
  L <- D - matrx
  # Compute Normalized Laplacian
  L <- as.matrix(L)
  D <- as.matrix(D)
  eig <- geigen(L,D, symmetric=T)
  rm(L,D)
  U <- eig$vectors[,1:rank]
  return(U)
}
