library(glmnet)
library(rTensor)

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


logi_deriv <- function(x, scale_par) {
  x[x <= -30] <- -30
  exp(-x / scale_par) / (scale_par * (1 + exp(-x / scale_par))^2)
}


proj_dist <- function(U, V) {
  norm(U %*% t(U) - V %*% t(V), "F")
}


grad_l_theta <- function(M1theta, M1A, type, scale_par) {
  eps <- 1e-15
  if (type == "logit") {
    temp <- logi_fun(M1theta, scale_par)
    return(logi_deriv(M1theta, scale_par) * ((1 - M1A) / (1 + eps - temp) - M1A / (eps + temp)))
  }
  if (type == "probit") {
    temp <- pnorm(M1theta / scale_par)
    return(dnorm(M1theta / scale_par) / scale_par * ((1 - M1A) / (1 - temp) - M1A / temp))
  }
  if (type == "poisson") {
    return(exp(M1theta) - M1A)
  }
}

obj_f <- function(M1A, M1theta, type, scale_par) {
  eps <- 1e-15
  if (type == "logit") {
    return(-sum(M1A * log(logi_fun(M1theta, scale_par) + eps) + (1 - M1A) * log(1 - logi_fun(M1theta, scale_par) + eps)))
  }
  if (type == "probit") {
    return(-sum(M1A * log(pnorm(M1theta / scale_par) + eps) + (1 - M1A) * log(1 - pnorm(M1theta / scale_par) + eps)))
  }
  if (type == "poisson") {
    return(sum(exp(M1theta) - M1A * M1theta))
  }
}

BIC_factor <- function(n, L, r, M) {
  return((2 * n * r - r * (r + 1) + L * M - 0.5 * M * (M + 1) + r * r * M) * log(n * n * L))
}

l_C_min <- function(U, V, W, M1A, cmax, type, rd, scale_par, sample_size) {
  r_1 <- ncol(U)
  Big_X <- kronecker(kronecker(W, V), U)
  vecA <- as.matrix(as.vector(M1A))
  if (rd == "rand") {
    ind_zero <- which(as.vector(M1A) == 0)
    ind_one <- setdiff(c(1:length(vecA)), ind_zero)
    ind <- c(sample(ind_zero, sample_size * 0.5), sample(ind_one, sample_size * 0.5))
    Big_X <- Big_X[ind, ]
    vecA <- vecA[ind, ]
  }
  # if (cmax==0){
  #   temp_df <- as.data.frame(cbind(vecA,Big_X/scale_par))
  #   result <- glm(V1~.-1, data=temp_df, family = binomial(link = type))
  #   return(M1C_min = matrix(result$coefficients, nrow = r_1)/10^6)
  # }else{
  result <- glmnet(Big_X / scale_par, vecA,
                   family = binomial(link = type),
                   alpha = 0, lambda = 0, intercept = FALSE,
                   lower.limits = -cmax, upper.limits = cmax
  )
  return(M1C_min = matrix((result$beta), nrow = r_1))
  # }
}


ProjectedGD <- function(Ini_list, Cmax, eta_outer,
                           tmax_outer, p_type, rd, show = TRUE, scale_par, sample_size)
{
  library(rTensor)
  library(glmnet)
  A <- Ini_list[[1]]
  U_int <- Ini_list[[2]]
  W_int <- Ini_list[[3]]
  delta <- Ini_list[[4]]
  data_modes <- A@modes
  r <- ncol(U_int)
  M <- ncol(W_int)
  core_modes <- c(r, r, M)
  M1A <- k_unfold(A, m = 1)@data
  ## Initialization for GD
  M1Ct <- l_C_min(U_int, U_int, W_int, M1A, Cmax, p_type, rd, scale_par, sample_size)
  M3Ct <- k_unfold(k_fold(M1Ct, m = 1, modes = core_modes), m = 3)@data
  ## Calculate Theta_0
  Up <- Ut <- U_int
  Wp <- Wt <- W_int
  M1thetat <- Ut %*% M1Ct %*% t(kronecker(Wt, Ut))
  obj <- obj_f(M1A, M1thetat, p_type, scale_par)
  obj_p <- obj
  cat("Running Grassmanian GD of", p_type, ":\n")
  for (t in 2:tmax_outer) {
    M1_l_grad <- grad_l_theta(M1thetat, M1A, p_type, scale_par)
    M3_l_grad <- k_unfold(k_fold(M1_l_grad, m = 1, modes = data_modes), m = 3)@data
    gradUt <- M1_l_grad %*% (kronecker(Wt, Ut) %*% t(M1Ct))
    gradWt <- M3_l_grad %*% (kronecker(Ut, Ut) %*% t(M3Ct))
    Ut <- svd(reg_fun(Ut - eta_outer * gradUt, delta[1]), nu = core_modes[1])[[2]]
    Wt <- svd(reg_fun(Wt - eta_outer * gradWt, delta[3]), nu = core_modes[3])[[2]]
    M1Ct <- l_C_min(Ut, Ut, Wt, M1A, Cmax, p_type, rd, scale_par, sample_size)
    M3Ct <- k_unfold(k_fold(M1Ct, m = 1, modes = core_modes), m = 3)@data
    M1thetat <- Ut %*% M1Ct %*% t(kronecker(Wt, Ut))
    if (show == TRUE) {
      cat("iteration", t, ": obj value:", obj, "\n")
      cat("iteration", t, ": UW_res:", UW_res = proj_dist(Ut, Up) + proj_dist(Wt, Wp), "\n")
    }
    obj <- obj_f(M1A, M1thetat, p_type, scale_par)
    if ((obj > obj_p) && (rd == "Non")) {
      break
    }
    obj_p <- obj
    Up <- Ut
    Wp <- Wt
  }
  if (rd == "rand") {
    M1Ct <- l_C_min(Ut, Ut, Wt, M1A, Cmax, p_type, "Non", scale_par, sample_size)
    M1thetat <- Ut %*% M1Ct %*% t(kronecker(Wt, Ut))
  }
  cat("BIC", 2 * obj + BIC_factor(data_modes[1], data_modes[3], r, M), "\n")
  list(
    "M1C" = M1Ct, "U_est" = Ut, "W_est" = Wt,
    "M1theta_est" = M1thetat
  )
}
