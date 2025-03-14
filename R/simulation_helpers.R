#  Functional Simulation: ---------
#' Simulate Dirichlet multinomial counts.
#'
#' @param Y0 Total sum of counts in DM
#' @param alpha
#'
#' @returns DM counts
#' @export
Dirichlet.multinomial <- function(Y0, alpha) {
  if (missing(Y0) || missing(alpha)) {
    stop("Total sum Y0 and/or alpha missing.")
  }

  # Create the data from the rmultinom
  dmData <- matrix(0, length(Y0), length(alpha))
  for (i in 1:length(Y0)) {
    # This is what I don't get... Why is it rmultinom + dirichlet?
    dmData[i, ] <- stats::rmultinom(1,
                                    Y0[i],
                                    dirmult::rdirichlet(1, alpha))
  }
  # Label the created data
  colnames(dmData) <- paste("Taxa", 1:ncol(dmData))

  return(dmData)
}




#' User defined correlation structure
#'
#' @param mi number of timepoints for ith sample
#' @param user_var 1: compound symmetry, 2: autoregressive, 3: indepdendent
#' @param structure Type of correlation 1 for (?) 2 for (?)
#' @param al ??
#'
#' @returns User defined correlation matrix.
#' @export
cor_user <- function(mi, user_var, structure, al) {
  if (structure == 1) {
    cor <- user_var * ((pracma::ones(mi, mi) - diag(rep(1, mi))) * al + diag(rep(1, mi)))
  }
  if (structure == 2) {
    need1 <- matrix(rep(1:mi, each = mi), nrow = mi)
    need2 <- matrix(rep(1:mi, each = mi), ncol = mi, byrow = TRUE)
    cor <- al^abs(need1 - need2) * user_var
  }
  if (structure == 3) {
    cor <- diag(mi) * user_var
  }
  return(cor)
}


#' Generating data simulation function
#'
#' @param n n constant of number of individuals
#' @param time vector of time values for each subject/time
#' @param miss_n missingness
#' @param cor_matrix correlation matrix
#' @param ranges ranges
#' @param alpha alphas for the Dirichlet multinomial
#' @param Z
#'
#' @returns Simulated data frame
#' @export
generate_data_longitudinal_compositional <- function(n,
                                                     time,
                                                     miss_n,
                                                     cor_matrix,
                                                     ranges,
                                                     alpha,
                                                     Z) {
  simulate_data_new <- list()
  for (i in 1:n) {
    simulate_data <- NULL
    # Set missingness
    miss_rate <- sample(miss_n, 1)
    newtime <- sort(sample(1:length(time), ceiling(miss_rate * length(time))))

    # Get the error matrix representing the longitudinal correlation
    error_matrix <- t(MASS::mvrnorm(n = ncol(alpha),
                                    mu = rep(0, length(time)),
                                    Sigma = cor_matrix,
                                    tol = 1e-6,
                                    empirical = FALSE,
                                    EISPACK = FALSE))

    for (j in 1:length(time)) {
      if (length(ranges) == 1) {
        Y0 <- ranges
      } else {
        Y0 <- sample(ranges, 1)
      }
      # Generate the compositional correlation
      dm <- Dirichlet.multinomial(Y0, as.numeric(alpha[j, ]))

      # Add the longitudinal correlation to the compositional correlation
      # Round it to be a count
      dm_witherror <- round(dm + error_matrix[j, ], 0)

      # Make sure no zeros.
      dm_witherror[dm_witherror < 0] <- 0

      # Define the
      total_n <- rowSums(dm_witherror)

      # Return data with the total n and id.
      dm_new <- c(dm_witherror, total_n, i)
      simulate_data <- rbind(simulate_data, dm_new)
    }

    # Include missingness.
    simulate_data_new[[i]] <- data.frame(time, simulate_data)[newtime, ]

    rownames(simulate_data_new[[i]]) <- NULL
  }
  return(simulate_data_new)
}

#' Generating data simulation function
#'
#' @param n n constant of number of individuals
#' @param time vector of time values for each subject/time
#' @param miss_n missingness
#' @param cor_matrix correlation matrix
#' @param ranges ranges
#' @param Z
#'
#' @returns Simulated data frame
#' @export
generate_data_longitudinal_compositionalZ <- function(n,
                                                      time,
                                                      miss_n,
                                                      cor_matrix,
                                                      ranges) {
  simulate_data_new <- list()
  for (i in 1:n) {
    simulate_data <- NULL

    # Three possible clusters with given
    fxn1 <- exp(cos(2 * pi * time))
    fxn2 <- exp(1 - 2 * exp(-6 * time))
    fxn3 <-  exp(-1.5 * time)

    Z <- stats::rbinom(n = length(time),
                       size = 1,
                       prob = 0.6)


    alpha <- data.frame(a1 = fxn1 - Z*fxn1,
                        a2 = fxn1 - Z*fxn1,
                        a3 = fxn1 - Z*fxn1,
                        a4 = fxn1 - Z*fxn1,
                        a5 = fxn2 - Z*(fxn2 - 1),
                        a6 = fxn2 - Z*(fxn2 - 1),
                        a7 = fxn2 - Z*(fxn2 - 1),
                        a8 = fxn2 - Z*(fxn2 - 1),
                        a9 =  fxn3 + Z*fxn1,
                        a10 = fxn3 + Z*fxn1,
                        a11 = fxn3 + Z*fxn1,
                        a12 = fxn3 + Z*fxn1)

    # Get the error matrix representing the longitudinal correlation
    error_matrix <- t(MASS::mvrnorm(n = ncol(alpha),
                                    mu = rep(0, length(time)),
                                    Sigma = cor_matrix,
                                    tol = 1e-6,
                                    empirical = FALSE,
                                    EISPACK = FALSE))

    for (j in 1:length(time)) {
      if (length(ranges) == 1) {
        Y0 <- ranges
      } else {
        Y0 <- sample(ranges, 1)
      }
      # Generate the compositional correlation
      dm <- Dirichlet.multinomial(Y0, as.numeric(alpha[j, ]))

      # Add the longitudinal correlation to the compositional correlation
      # Round it to be a count
      dm_witherror <- round(dm + error_matrix[j, ], 0)

      # Make sure no zeros.
      dm_witherror[dm_witherror < 0] <- 0

      # Define the
      total_n <- rowSums(dm_witherror)

      # Return data with the total n and id.
      dm_new <- c(dm_witherror, total_n, i, Z[j])
      simulate_data <- rbind(simulate_data, dm_new)
    }

    # Include missingness.
    # Set missingness
    miss_rate <- sample(miss_n, 1)
    newtime <- sort(sample(1:length(time), ceiling(miss_rate * length(time))))
    simulate_data_new[[i]] <- data.frame(time, simulate_data)[newtime, ]

    rownames(simulate_data_new[[i]]) <- NULL
  }
  return(simulate_data_new)
}

#' Use Chenyangs setup to simulate count data wtih 3 clusters
#'
#' No external variables.
#'
#' @returns data Matrix with columns time, individual, capture number, totaln, counts
#' @export
sim_noZ <- function(n = 20,
                    range_start = 5000,
                    range_end = 20000,
                    nknots = 3,
                    order = 3,
                    user_var = 1000,
                    structure = 1,
                    al = 0.4
                    ){
  # Time points are a sequence between 0 and 1
  time <- seq(0, 1, 0.05)

  # Three possible clusters with given
  fxn1 <- exp(cos(2 * pi * time))
  fxn2 <- exp(1 - 2 * exp(-6 * time))
  fxn3 <-  exp(-1.5 * time)

  curve_paramenter <- data.frame(a1 = fxn1,
                                a2 = fxn1,
                                a3 = fxn1,
                                a4 = fxn1,
                                a5 = fxn2,
                                a6 = fxn2,
                                a7 = fxn2,
                                a8 = fxn2,
                                a9 = fxn3,
                                a10 = fxn3,
                                a11 = fxn3,
                                a12 = fxn3)

  # Simulation parameters.
  ranges <- range_start:range_end
  K <- ncol(curve_paramenter)

  # How many missing times?
  miss_n <- seq(0.6, 1, by = 0.1)

  # Set up the longitudinal correlation structure
  # What are these numbers?
  cor_matrix <- cor_user(length(time),
                         user_var = user_var,
                         structure = structure,
                         al = al)

  # Simulate data.
  generate.data <- generate_data_longitudinal_compositional(n,
                                                            time,
                                                            miss_n,
                                                            cor_matrix,
                                                            ranges,
                                                            curve_paramenter)
  # Organize data.
  data <- dplyr::bind_rows(generate.data)
  names(data) <- c("time",
                   paste0("Taxa.", 1:K),
                   "total_n",
                   "individual")

  data$Capture.Number <- data$time * (length(unique(data$time)) - 1) + 1
  data$individual <- as.factor(data$individual)
  return(data)
}



#' Use Chenyangs setup to simulate count data wtih 3 clusters
#'
#' No external variables.
#'
#' @returns data Matrix with columns time, individual, capture number, totaln, counts
#' @export
sim_Z_longitudinal <- function(n = 20,
                               range_start = 5000,
                               range_end = 20000,
                               nknots = 3,
                               K = 12,
                               order = 3,
                               user_var = 1000,
                               structure = 3,
                               al = 0.4,
                               miss_p = 0.6
){
  # Time points are a sequence between 0 and 1
  time <- seq(0, 1, 0.05)

  # Simulation parameters for setting Y0
  ranges <- range_start:range_end

  # How many missing times?
  miss_n <- seq(miss_p, 1, by = 0.1)

  # Set up the longitudinal correlation structure
  # What are these numbers?
  cor_matrix <- cor_user(length(time),
                         user_var = user_var,
                         structure = structure,
                         al = al)

  # Simulate data.
  generate.data <- generate_data_longitudinal_compositionalZ(n,
                                                            time,
                                                            miss_n,
                                                            cor_matrix,
                                                            ranges)
  # Organize data.
  data <- dplyr::bind_rows(generate.data)
  names(data) <- c("time",
                   paste0("Taxa.", 1:K),
                   "total_n",
                   "individual",
                   "Z")

  data$Capture.Number <- data$time * (length(unique(data$time)) - 1) + 1
  data$individual <- as.factor(data$individual)

  return(data)
}




# Simulate from B-spline --------------------------------------------------


#' Simulate timepoints
#'
#' Get a collection of timepoints along the same range
#'
#' @param n number of subjects
#' @param maxt maximum time
#' @param max_mi maximum number of timepoints for each subject
#'
#' @returns list
#' @export
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
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Formatted Z matrix
#' @export
sim_Z <- function(mi_vec){

  M <- sum(mi_vec)
  Z0 <- rep(1, M)
  Z1 <- stats::rnorm(M)
  Z2 <- stats::rbinom(M, 1, 0.5)
  Z <- matrix(c(Z0, Z1, Z2), ncol = 3)
  return(Z)
}


#' Simulate Yij counts
#'
#' Use Dirichlet Multinomial distribution to simulate counts
#'
#' @param i subject index
#' @param j time index
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param Y_ij0 Total sum of counts across all K
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Y_i
#' @export
sim_Y_ij <- function(i, j, beta, Z, B, Y_ij0, K, mi_vec, i_index){
  Y_ij <- MGLM::rdirmn(n = 1,
                       size = Y_ij0,
                       alpha = get_alpha_ij(i = i,
                                            j = j,
                                            beta = beta,
                                            Z = Z,
                                            B = B,
                                            K = K,
                                            i_index = i_index))
  return(Y_ij) # Returns a 1xK matrix
}

# Simulate Y for all i, j, k
sim_Yi <- function(i, beta, Z, B, K, mi_vec, i_index){
  mi <- mi_vec[i]
  Y_i <- matrix(nrow = mi, ncol = K)
  for (j in 1:mi) {
    Y_ij0 <- sample(100:500, 1)
    Y_i[j,] <- sim_Y_ij(i, j,
                        beta, Z, B,
                        Y_ij0, K,
                        mi_vec, i_index)
  }
  return(Y_i)
}


#' Simulate all Y counts
#'
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param B B spline basis matrix of dimension (N x P)
#' @param K Number of responses
#' @param mi_vec vector of the number of timepoints for each sample. Of length n
#'
#' @returns Y
#' @export
sim_Y <- function(beta, Z, B, K, mi_vec, i_index){
  Y <- matrix(nrow = 0, ncol = K)
  n <- length(mi_vec)

  for (i in 1:n) {
    Y <- rbind(Y,
               sim_Yi(i,
                      beta,
                      Z,
                      B,
                      K,
                      mi_vec,
                      i_index))
  }
  return(Y)
}


#' Title
#'
#' @param seed Seed for the random number generator
#'
#' @returns list
#' @export
base_sim <- function(seed = 124){
  set.seed(seed)
  time_list <- sim_timepoints(n = 5)
  time <- time_list$X$time
  mi_vec <- time_list$mi_vec
  i_index <- c(0, cumsum(mi_vec))
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
             mi_vec = mi_vec,
             i_index = i_index)

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
