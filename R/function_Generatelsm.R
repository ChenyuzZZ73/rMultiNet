library(rTensor)
library(glmnet)
library(truncnorm)
library(lpSolve)
library(matlabr)
#####################################
logi_fun <- function(x, scale_par) {
  1 / (1 + exp(-x / scale_par))
}

reg_fun <- function(U, delta) {
  U <- as.matrix(U)
  # print(U)
  ind <- sqrt(rowSums(U^2)) > delta
  if (sum(ind) > 1) {
    U[ind, ] <- U[ind, ] * delta / sqrt(rowSums(U[ind, ]^2))
  } else if (sum(ind) == 1) {
    U[ind, ] <- U[ind, ] * delta / sqrt(sum(U[ind, ]^2))
  }
  return(U)
}


#############################################
generate_U <- function(Umean, n, d) {
  U <- matrix(0, nrow = n, ncol = d)
  for (i in 1:n) {
    U[i, ] <- rnorm(d, mean = Umean, sd = 1)
  }
  U <- svd(U)$u
  U <- sqrt(n) * U
  return(U)
}

generate_C <- function(cmax, d, M, int_type) {
  if(int_type == "Uniform")
  {
    C <- array(runif(d * d * M, -cmax, cmax), dim = c(d, d, M))
    for (i in 1:M) {
      A <- matrix(runif(d * d, -cmax, cmax), nrow = d)
      C[, , i] <- A %*% t(A)
    }
  }
  if (int_type == "Norm")
  {
    C <- array(rnorm(d*d*M,0,cmax), dim = c(d,d,M))
    for (i in 1:M)
    {
      A <- matrix(rnorm(d*d*M,0,cmax), nrow = d)
      C[, , i] <- A %*% t(A)
    }
  }
  return(C)
}

#############################################
generate_lsm <- function(U, C, type, scale_par) {
  Theta <- U %*% C %*% t(U)
  if (type == "logit")
  {
    P <- logi_fun(Theta, scale_par)
  }
  if (type == "probit")
  {
    P <- pnorm(Theta / scale_par)
  }
  n <- nrow(U)
  A <- matrix(NA, ncol = n, nrow = n)
  diag(A) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n)
    {
      A[i, j] <- rbinom(1, 1, P[i, j])
      A[j, i] <- A[i, j]
    }
  }
  return(list(A, Theta))
}


Generatelsm <- function(n, L, M, d, U_mean, cmax, int_type, kernel_fun, scale_par)
{
  Pi_network_type=rep(1/M,M);
  ## Generate network labels
  network_type=sample(x = M, size = L, replace = TRUE, prob = Pi_network_type)
  while (length(unique(network_type))!=M) {
    network_type=sample(x = M, size = L, replace = TRUE, prob = Pi_network_type)
  }
  my_tensor <- array(0, dim = c(n, n, L))
  # my_U <- array(0, dim = c(n,d,M))
  my_Theta <- array(0, dim = c(n, n, L))
  my_W <- matrix(0, nrow = L, ncol = M)
  U <- generate_U(U_mean,n,d)
  C <- generate_C(cmax,d,M, int_type)
  for (l in 1:L) {
    gen_list <- generate_lsm(U, C[, , network_type[l]], kernel_fun, scale_par)
    my_tensor[, , l] <- gen_list[[1]]
    library(rTensor)
    adj_tensor <- as.tensor(my_tensor)
    # my_U[,,network_type[l]]=gen_list[[2]]
    my_Theta[, , l] <- gen_list[[2]]
    Theta <- as.tensor(my_Theta)
    my_W[l, network_type[l]] <- 1
  }
  return(list(adj_tensor, Theta , my_W))
}
