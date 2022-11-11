# Given n by n similarity this function first calculate the Laplacian
# matrix L then generate n by ncol matrix U of top ncol eigenvectors of L.
#
# Args:
#     matrx: an n by n matrix
#
#     ncol: number of columns of the output matrix U
#
# Returns:
#    U: n by ncol numeric matrix that contains the ncol tops
#       eigenvectors of Laplacian matrix as column
#
SpecClustering <- function(tnsr, ncol, embedding_type = "Layer")
{
  library(Matrix)
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
  U <- eig$vectors[,1:ncol]
  return(U)
}
