Community_cluster_km <- function(decomp,type,cluster_number)
{
    if(type == "N")
   {
     embedding = decomp[["U"]][[3]]
     re = kmeans(embedding[,c(2,3)],cluster_number,iter.max = 100)
     png("cluster_network.png")
     plot(embedding[,c(2,3)], col=as.factor(re$cluster),main = "cluster result", xlab = "second eigen vector", ylab = "third eigen vector")
     dev.off()
    }
    if(type == "n")
    {
      embedding = decomp[["U"]][[1]]
      #library(igraph)
      library(plotly)
      re = kmeans(embedding[,c(2,3,4)],cluster_number, iter.max = 100)
      plot_ly(x = embedding[,2], y =embedding[,3], z =embedding[,4],size=2,color=as.character(re$cluster))
    }

}


Community_cluster_dbscan <- function(decomp,type,eps_value,pts_value)
{
  if(type == "N")
  {
    embedding = decomp[["U"]][[3]]
    dbscan_tensor=dbscan(embedding[,2:3],eps=eps_value, minPts =pts_value)
    png("cluster_network.png")
    plot(embedding[,2:3], col=as.factor(dbscan_tensor$cluster),main = "cluster result", xlab = "second eigen vector", ylab = "third eigen vector")
    dev.off()
  }
  if(type == "n")
  {
    embedding = decomp[["U"]][[1]]
    #library(igraph)
    library(plotly)
    dbscan_tensor=dbscan(embedding[,c(2,3,4)],eps=eps_value,minPts = eps_value)
    plot_ly( x = embedding[,2], y =embedding[,3], z =embedding[,4],size=2
             ,color=as.character(dbscan_tensor$cluster))

  }

}
