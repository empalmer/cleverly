#' Simulate timepoints
#'
#' Get a collection of timepoints along the same range
#'
#' @param n
#' @param maxt
#' @param max_mi
#'
#' @returns
#' @export
#'
#' @examples
sim_timepoints <- function(n, maxt = 9, max_mi = 4){
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


#' Simulate 2 external variables, one continuous and one binary
#'
#' @param mi_vec
#'
#' @returns
#' @export
#'
#' @examples
sim_Z <- function(mi_vec){

  M <- sum(mi_vec)
  Z0 <- rep(1, M)
  Z1 <- rnorm(M)
  Z2 <- rbinom(M, 1, 0.5)
  Z <- matrix(c(Z0, Z1, Z2), ncol = 3)
  return(Z)
}


#' Simulate Yij counts
#'
#' Use Dirichlet Multinomial distribution to simulate counts
#'
#' @param i
#' @param j
#' @param beta
#' @param Z
#' @param B
#' @param Y_ij0
#' @param K
#' @param mi_vec
#'
#' @returns
#' @export
#'
#' @examples
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


#' Simulate all Y counts
#'
#' @param beta
#' @param Z
#' @param B
#' @param K
#' @param mi_vec
#'
#' @returns
#' @export
#'
#' @examples
sim_Y <- function(beta, Z, B, K, mi_vec){
  Y <- matrix(nrow = 0, ncol = K)
  n <- length(mi_vec)
  for (i in 1:n) {
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

#' Simulation from B-spline
#'
#' simulate counts and external variables directly from b-spline coefficients
#'
#' @param seed
#'
#' @returns
#' @export
#'
#' @examples
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

# Chenyang's data simulation functions: ---------
Dirichlet.multinomial <- function(Nrs, shape) {
  if (missing(Nrs) || missing(shape)) {
    stop("Nrs and/or shape missing.")
  }

  # Create the data from the rmultinom
  dmData <- matrix(0, length(Nrs), length(shape))
  for (i in 1:length(Nrs)) {
    dmData[i, ] <- stats::rmultinom(1, Nrs[i], dirmult::rdirichlet(1, shape))
  }

  # Label the created data
  colnames(dmData) <- paste("Taxa", 1:ncol(dmData))

  return(dmData)
}


rep.row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}

rep.col <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}


combine_data <- function(generate.data, indiv_n) {
  simulate_data <- NULL
  for (i in 1:indiv_n) {
    simulate_data <- rbind(simulate_data, generate.data[[i]])
  }
  return(simulate_data)
}

cor_user <- function(mi, user_var, structure, al) {
  if (structure == 1) {
    return(user_var * ((pracma::ones(mi, mi) - diag(rep(1, mi))) * al + diag(rep(1, mi))))
  }
  if (structure == 2) {
    return(al^abs(rep.col(1:mi, mi) - rep.row(1:mi, mi)) * user_var)
  }
}


generate_data_new <- function(indiv_n, time, miss_n, cor_matrix, ranges, curveparamenter) {
  simulate_data_new <- list()
  for (i in 1:indiv_n) {
    Get_alpha <- curveparamenter
    simulate_data <- NULL
    var_yij_final <- NULL
    miss_rate <- sample(miss_n, 1)
    newtime <- sort(sample(1:length(time), ceiling(miss_rate * length(time))))
    error_matrix <- t(MASS::mvrnorm(n = ncol(curveparamenter),
                                    mu = rep(0, length(time)),
                                    Sigma = cor_matrix,
                                    tol = 1e-6,
                                    empirical = FALSE,
                                    EISPACK = FALSE))

    for (j in 1:length(time)) {
      if (length(ranges) == 1) {
        total_n_old <- ranges
      } else {
        total_n_old <- sample(ranges, 1)
      }
      dm <- Dirichlet.multinomial(total_n_old, as.numeric(Get_alpha[j, ]))
      dm_witherror <- round(dm + error_matrix[j, ], 0)
      dm_witherror[dm_witherror < 0] <- 0
      total_n <- rowSums(dm_witherror)
      dm_new <- c(dm_witherror, total_n, i)
      simulate_data <- rbind(simulate_data, dm_new)
    }


    simulate_data_new[[i]] <- data.frame(time, simulate_data)[newtime, ]

    rownames(simulate_data_new[[i]]) <- NULL
  }
  return(simulate_data_new)
}



combine_data <- function(generate.data, indiv_n){
  simulate_data <- NULL
  for (i in 1:indiv_n) {
    simulate_data <- rbind(simulate_data,generate.data[[i]])
  }
  return(simulate_data)
}

#' Use Chenyangs setup to simulate count data wtih 3 clusters
#'
#' No external variables.
#'
#' @returns
#' @export
#'
#' @examples
sim_noZ <- function(){
  time <- seq(0,1,0.05)
  indiv_n <- 20
  a1 <- exp(cos(2 * pi * time))
  a2 <- exp(cos(2 * pi * time))
  a3 <- exp(cos(2 * pi * time))
  a4 <- exp(cos(2 * pi * time))

  a5 <- exp(1 - 2 * exp(-6 * time))
  a6 <- exp(1 - 2 * exp(-6 * time))
  a7 <- exp(1 - 2 * exp(-6 * time))
  a8 <- exp(1 - 2 * exp(-6 * time))

  a9  <- exp(-1.5 * time)
  a10 <- exp(-1.5 * time)
  a11 <- exp(-1.5 * time)
  a12 <- exp(-1.5 * time)

  curveparamenter <- data.frame(a1, a2,  a3,  a4,
                                a5, a6,  a7,  a8,
                                a9, a10, a11, a12)

  ranges <- 5000:20000
  nknots <- 3
  order <- 3
  n <- 12
  gamma1 <- 1/300

  miss_n <- seq(0.6, 1, by = 0.1)
  cor_matrix <- cor_user(length(time),
                         user_var = 1000,
                         structure = 1,
                         al = 0.4)
  generate.data <- generate_data_new(indiv_n,
                                     time,
                                     miss_n,
                                     cor_matrix,
                                     ranges,
                                     curveparamenter)

  simulate_data_new <- combine_data(generate.data, indiv_n)
  names(simulate_data_new) <- c("time",
                                "Taxa.1","Taxa.2",
                                "Taxa.3","Taxa.4",
                                "Taxa.5","Taxa.6",
                                "Taxa.7","Taxa.8",
                                "Taxa.9","Taxa.10",
                                "Taxa.11","Taxa.12",
                                "total_n","individual")
  simulate_data_new$Capture.Number <- simulate_data_new$time * (length(unique(simulate_data_new$time)) - 1) + 1
  simulate_data_new$individual <- as.factor(simulate_data_new$individual)
  data <- simulate_data_new

  return(data)
}
