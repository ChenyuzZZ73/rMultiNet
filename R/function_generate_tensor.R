###Function generate_tesnor: generates MMSBM
###Input:
###      n: number of nodes
###      m: types of networks
###      K: number of communities
###      L: number of layers
###      d: average degree (a L-dimensional vector)
###      r: out-in ratio (q/p, a L-dimensional vector)
###      Pi: community size (defult as 1/K, a K-dimensional vector)
###      Pi_network_type: size of network types (defult as 1/m, an m-dimensional vector)
###Output: list(arrT(recording all the neworks),global membership,local membership,network type)

require(rTensor)

generate_sbm<-function(n,K,r,d,membership)
{
  p=d*K/(n*(1+(K-1)*r))
  q=p*r
  B=matrix(q,K,K)
  diag(B)=p

  Z=matrix(0,nrow=n,ncol=K)
  for(i in 1:n)
  {
    Z[i,membership[i]]=1
  }

  P=Z%*%B%*%t(Z)

  A=matrix(NA,ncol=n,nrow=n)
  diag(A)=0
  for (i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      A[i,j]=rbinom(1,1,P[i,j])
      A[j,i]=A[i,j]
    }
  }
  return(A)
}

generate_tensor<-function(n,m,K,L,r,d,Pi,Pi_network_type)
{
  ###Generate Membership matrix n by m
  local_membership=matrix(NA,nrow=n,ncol=m)
  colnames(local_membership)=c(1:m)
  for(i in 1:m)
  {
    local_membership[,i]=sample(x = K, size = n, replace = TRUE, prob = Pi)
  }

  ###Assign network type for each layer
  network_type=sample(x = m, size = L, replace = TRUE, prob = Pi_network_type)

  mylist=list()

  for (i in 1:L)
  {
    mylist[[i]]=generate_sbm(n=n,K=K,r=r[i],d=d[i],membership=local_membership[,network_type[i]])
  }

  global_membership=rep(0,n)
  for(i in 1:m)
  {
    global_membership=global_membership+K^(i-1)*(local_membership[,i]-1)
  }
  indices <- c(n,n,L)
  arr<-array(as.numeric(unlist(mylist)), dim=indices)
  arrT <- as.tensor(arr)
  output<-list(arrT,global_membership,local_membership,network_type)
  return(output)
}
