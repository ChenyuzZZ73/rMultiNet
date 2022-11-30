#' @title      A function for embedding visuilization
#' @name       Embedding_network
#' @description a method to visuilize the result of embedding
#' @param      network_membership
#' @param      L  the number of network membership(n or L)
#' @param      paxis = 2    the  number of axis to plot
#'@import      ggplot2
#' @export

Embedding_network <- function(network_membership,L, paxis = 2)
{
  included = rep(NA,L)
  for(i in 1:L)
  {
    included[i] = i
  }
  #network_membership = decomp[["U"]][[3]]
  d_2 = dim(network_membership)
  if(paxis>2)
  {
    data <- as.data.frame(network_membership[,2:paxis])
    return(data)
  }
  else{
    png(file = "embedding_network.png")
    plot(x = network_membership[,2],y = network_membership[,3], xlab="second eigen vector", ylab="third eigen vector")
    text(network_membership[,2],network_membership[,3],as.character(included),cex = 1,col = "red")
    dev.off()
  }
}

