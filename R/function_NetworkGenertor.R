source(file = "R/function_generate_tensor.R")

MMSBM_generation <- function(n,m,L,K,d,r)
{

  ###Initialize the average degree list
  if (is.null(d))
  {
    d_list = abs(rnorm(L,5,5))
  }
  else
    d_list = rep(d,L)

  ###Initialize the out-in ratio of each layer list
  if (is.null(r))
  {
    r_list = rep(0.4,L)
  }
  else
    r_list = rep(r,L)
  Pi_network_type = rep(1/m,m)
  Pi_network_community = rep(1/K,K)
  library(rTensor)
  temp = generate_tensor(n=n,m=m,K=K,L=L,r=r_list,d=d_list,Pi=Pi_network_community,Pi_network_type=Pi_network_type)
  arrT = temp[[1]]
  Global_node_type_power = temp[[2]]
  number_gloab_community = length(unique(Global_node_type_power))
  return(arrT)

}
