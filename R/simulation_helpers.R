

sim_timepoints <- function(n, maxt = 9, max_mi = 4){
  # generate timepoints: ----
  mis <- sample(3:max_mi, n, replace = TRUE)

  sample_ids <- numeric()
  timepoints <- numeric()

  for (i in 1:n) {
    timepoints <- c(timepoints, sort(sample(1:maxt, mis[i])))
    sample_ids <- c(sample_ids, rep(paste0("ID", i), mis[i]))
  }
  return(list(X = data.frame(id = sample_ids, time = timepoints),
              mis = mis))
}



sim_Z <- function(mis){
  # External variables: ------
  M <- sum(mis)


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
sim_Y_ij <- function(i, j, beta, Z, B, Y_ij0, K, mi){
  Y_ij <- MGLM::rdirmn(n = 1,
                       size = Y_ij0,
                       alpha = get_alpha_ij(i = i,
                                            j = j,
                                            beta = beta,
                                            Z = Z,
                                            B = B,
                                            K = K,
                                            mi = mi))
  return(Y_ij) # Returns a 1xK matrix
}

# Simulate Y for all i, j, k
sim_Yi <- function(i, beta, Z, B, K, mi){
  Y_i <- matrix(nrow = mi, ncol = K)
  for (j in 1:mi) {
    Y_ij0 <- sample(50:150, 1)
    Y_i[j,] <- sim_Y_ij(i, j, beta, Z, B, Y_ij0, K, mi)
  }
  return(Y_i)
}


sim_Y <- function(beta, Z, B, K, mis){
  #Y <- matrix(nrow = sum(mis), ncol = K)
  Y <- matrix(nrow = 0, ncol = K)
  n <- length(mis)
  for (i in 1:n) {
    #Y[((i - 1)*mis[i] + 1):(i*mis[i]), ] <- sim_Yi(i, beta, Z, B)
    Y <- rbind(Y, sim_Yi(i, beta, Z, B, K, mis[i]))
  }
  return(Y)
}

