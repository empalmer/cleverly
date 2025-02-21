

sim_timepoints <- function(n, maxt = 9, max_mi = 4){
  # generate timepoints: ----
  mi_vec <- sample(3:max_mi, n, replace = TRUE)

  sample_ids <- numeric()
  timepoints <- numeric()

  for (i in 1:n) {
    timepoints <- c(timepoints, sort(sample(1:maxt, mi_vec[i])))
    sample_ids <- c(sample_ids, rep(paste0("ID", i), mi_vec[i]))
  }
  return(list(X = data.frame(id = sample_ids, time = timepoints),
              mi_vec = mi_vec))
}



sim_Z <- function(mi_vec){
  # External variables: ------
  M <- sum(mi_vec)


  Z0 <- rep(1, M)
  Z1 <- rnorm(M)
  Z2 <- rbinom(M, 1, 0.5)

  #Z_df <- data.frame(id = sample_ids, time = timepoints_i, Z0, Z1, Z2)
  Z <- matrix(c(Z0, Z1, Z2), ncol = 3)

  return(Z)
}



# # Set beta
# # For cluster 1 (k = 1, 2)
# betaC1 <- matrix(c(rep(c(1, 1, 1, 1, 1, 1), 2),     #l = 0
#                    rep(c(-2, -2, -2, -2, -2, -2),2),  #l = 1
#                    rep(c(1, 2, 3, 4, 5, 6), 2)), ncol = 3)  #l = 2
#
# betaC2 <- -betaC1
#
# beta <- rbind(betaC1, betaC2)
#
#
#
# # Simulate Y ----------------------------------------------------
#
#
sim_Y_ij <- function(i, j, beta, Z, B, Y_ij0, K, mi_vec){
  Y_ij <- MGLM::rdirmn(n = 1,
                       size = Y_ij0,
                       alpha = get_alpha_ij(i = i,
                                            j = j,
                                            beta = beta,
                                            Z = Z,
                                            B = B,
                                            K = K,
                                            mi_vec = mi_vec))
  return(Y_ij) # Returns a 1xK matrix
}

# Simulate Y for all i, j, k
sim_Yi <- function(i, beta, Z, B, K, mi_vec){
  mi <- mi_vec[i]
  Y_i <- matrix(nrow = mi, ncol = K)
  for (j in 1:mi) {
    Y_ij0 <- sample(50:150, 1)
    Y_i[j,] <- sim_Y_ij(i, j, beta, Z, B, Y_ij0, K, mi_vec)
  }
  return(Y_i)
}


sim_Y <- function(beta, Z, B, K, mi_vec){
  #Y <- matrix(nrow = sum(mi_vec), ncol = K)
  Y <- matrix(nrow = 0, ncol = K)
  n <- length(mi_vec)
  for (i in 1:n) {
    #Y[((i - 1)*mi_vec[i] + 1):(i*mi_vec[i]), ] <- sim_Yi(i, beta, Z, B)
    Y <- rbind(Y,
               sim_Yi(i,
                      beta,
                      Z,
                      B,
                      K,
                      mi_vec))
  }
  return(Y)
}

base_sim <- function(seed = 124){
  set.seed(seed)
  time_list <- sim_timepoints(n = 5)
  time <- time_list$X$time
  mi_vec <- time_list$mi_vec
  is <- rep(1:5, mi_vec)
  M <- sum(mi_vec)
  K <- 4
  Z <- sim_Z(mi_vec)
  B <- get_B(time, order = 3, nknots = 3)
  P <- 6

  # 2 clusters
  betaC1 <- matrix(c(rep(c(1, 1, 1, 1, 1, 1), 2),     #l = 0
                     rep(c(-2, -2, -2, -2, -2, -2),2),  #l = 1
                     rep(c(1, 2, 3, 4, 5, 6), 2)), ncol = 3)#l = 2
  betaC2 <- .5*betaC1
  beta <- rbind(betaC1, betaC2)

  Y <- sim_Y(beta = beta,
             Z = Z,
             B = B,
             K = K,
             mi_vec = mi_vec)

  return(list(Y = Y,
              Z = Z,
              B = B,
              is = is,
              beta = beta,
              mi_vec = mi_vec,
              time = time,
              K = K,
              P = P))
}
