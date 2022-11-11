###Function InitializationMMSBM: initialization for HOSVD
###Input:
###      tnsr

###Output:list of U1, U2 and U3

library(rTensor)

InitializationMMSBM<-function(tnsr, ranks=NULL)
{
  num_modes <- tnsr@num_modes
  U_list <- vector("list", num_modes)
  temp_mat <-matrix(0,ncol=tnsr@modes[1], nrow=tnsr@modes[2])
  for(i in 1:tnsr@modes[3])
  {
    temp_mat=temp_mat+tnsr@data[,,i]
  }
  U_list[[1]] <- eigen(temp_mat,symmetric = T)$vector[,c(1:ranks[1])]
  U_list[[2]] <- U_list[[1]]
  outer_production=NULL
  for(i in 1:dim(U_list[[1]])[1])
  {
    row_matrix=NULL
    for (j in 1:dim(U_list[[1]])[2])
    {
      temp=U_list[[1]][i,j]*U_list[[2]]
      row_matrix=cbind(row_matrix,temp)
    }
    outer_production=rbind(outer_production,row_matrix)
  }
  temp_mat <- rs_unfold(tnsr, m = 3)@data %*% outer_production
  U_list[[3]] <- svd(temp_mat, nu = ranks[3])$u
  return(U_list)
}
