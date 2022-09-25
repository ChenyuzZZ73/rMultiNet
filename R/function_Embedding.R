library(rTensor)
source(file = "R/function_Initialization.R")
source(file = "R/function_PowerIteration.R")

Embedding_network <- function(decomp,L)
{
  included = rep(NA,L)
  for(i in 1:L)
  {
    included[i] = i
  }
  network_membership = decomp[["U"]][[3]]
  core_tensor = decomp[["Z"]]

  d_2 = dim(network_membership)
  png(file = "embedding_network.png")
  plot(network_membership[,2:3],xlab="second eigen vector", ylab="third eigen vector")
  text(network_membership[,2:3],as.character(included),cex = 1,col = "red")
  dev.off()
}

Embedding_nodes <- function(decomp,N)
{
  included = rep(NA,N)
  for(i in 1:N)
  {
    included[i] = i
  }

  nodes_membership = decomp[["U"]][[1]]
  core_tensor = decomp[["Z"]]
  d_1 = dim(nodes_membership)
  png(file="embedding_nodes.png")
  plot(nodes_membership[,2:3],xlab="second eigen vector", ylab="third eigen vector")
  text(nodes_membership[,2:3],as.character(included),cex = 1,col = "green")
  dev.off()
}
