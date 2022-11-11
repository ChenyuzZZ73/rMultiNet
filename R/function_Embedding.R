require(rTensor)
require(ggplot2)

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
  }
  else{
    png(file = "embedding_network.png")
    plot(network_membership[,2:3],xlab="second eigen vector", ylab="third eigen vector")
    text(network_membership[,2:3],as.character(included),cex = 1,col = "red")
    dev.off()
  }
}

