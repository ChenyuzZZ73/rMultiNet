library(rTensor)
## Intialization for U,W
Init_UW <- function(tensor_A, n, r, M, Utrue, Wtrue, pertub = 0.1, int_type) {
  L <- tensor_A@modes[3]
  M1A <- k_unfold(tensor_A, m = 1)@data
  M3A <- k_unfold(tensor_A, m = 3)@data
  if(int_type == "spec")
  {
    U_int <- svd(M1A, nu = r)[["u"]]
    W_int <- svd(M3A, nu = M)[["u"]]
    return(list("U_int" = U_int, "W_int" = W_int))
  }
  if (int_type == "rand") {
    # random Initialization
    U_int <- generate_U(0, n, r) / sqrt(n)
    W_int <- svd(M3A, nu = M)[["u"]]
    return(list("U_int" = U_int, "W_int" = W_int))
  }
  if (int_type == "warm") {
    # warm Initialization
    per_U <- matrix(runif(n * n, -0.1, 0.1), ncol = n)
    per_W <- matrix(runif(L * L, -pertub, pertub), ncol = L)
    U_int <- svd(Utrue %*% t(Utrue) / n + per_U, nu = r)[["u"]]
    W_int <- svd(Wtrue %*% t(Wtrue) + per_W, nu = M)[["u"]]
    return(list("U_int" = U_int, "W_int" = W_int))
  }
}


Initilization <- function(gen_list, n, r, M, pertub = 0.1, int_type)
{
   adj_tensor <- gen_list[[1]]
   theta <- gen_list[[2]]
   tuc_res <- tucker(theta, ranks = c(r,r,M), max_iter = 25, tol = 1e-05)
   Utrue_scaled <- tuc_res$U[[1]]
   Wtrue <- tuc_res$U[[3]]
   delta_1 <- sqrt(max(rowSums(Utrue_scaled^2)))
   delta_3 <- sqrt(max(rowSums(Wtrue^2)))
   delta <- c(delta_1,delta_1,delta_3)
   initial_UW <- Init_UW(adj_tensor,n,r,M, Utrue_scaled, Wtrue, pertub, int_type)
   U_int <- initial_UW$U_int;
   W_int <- initial_UW$W_int;
   return(list(adj_tensor,U_int, W_int, delta))
}
