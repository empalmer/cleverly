
# Simulate the data (DM) --------------------------------------------------

#' Simulate data
#'
#' Use Chenyang's setup to simulate count data with 3 clusters,
#' one external variable binary with prob prob1 of being 1
#'
#' @param n Number of subjects
#' @param range_start Starting count range for the Dirichlet multinomial
#' @param range_end Ending count range for the Dirichlet multinomial
#' @param nknots Number of knots
#' @param K Number of responses
#' @param order Order of the B-spline basis
#' @param user_var User supplied variance
#' @param cor_str Type of correlation structure (IND, CON, AR1)
#' @param rho Correlation parameter
#' @param miss_p proportion of missing samples
#' @param prob1 Equivalent to p if Z ~ Binomial(p = prob1)
#' @param Z_type Type of external variable (binary or continuous)
#' @param maxt Time scale goes from 0 to maxt
#' @param baseline_fxns list of functions to use for the baseline
#' @param slope_fxns list of functions to use for the slope
#'
#' @returns data Matrix with columns time, individual, capture number, totaln, counts
#' @export
simulation_data <- function(n = 20,
                            K = 12,
                            range_start = 5000,
                            range_end = 20000,
                            Z_type = "binary",
                            user_var = 1000,
                            cor_str,
                            nknots = 3,
                            order = 3,
                            rho = 0.4,
                            miss_p = 0.6,
                            prob1 = 0.4,
                            maxt = 1,
                            baseline_fxns = list(
                              function(t) 2 * cos(2 * pi * t),
                              function(t) 2 * cos(2 * pi * t),
                              function(t) 2 * cos(2 * pi * t),
                              function(t) 2 * cos(2 * pi * t),
                              function(t) 1 - 2 * exp(-6 * t),
                              function(t) 1 - 2 * exp(-6 * t),
                              function(t) 1 - 2 * exp(-6 * t),
                              function(t) 1 - 2 * exp(-6 * t),
                              function(t) -3 * t + 3,
                              function(t) -3 * t + 3,
                              function(t) -3 * t + 3,
                              function(t) -3 * t + 3
                            ),
                            # Slope functions
                            slope_fxns = list(
                              function(t) sqrt(t),
                              function(t) 2 * t,
                              function(t) -3 * t + 3,
                              function(t) 4 * t,
                              function(t) -5 * t + 5,
                              function(t) -2 * t,
                              function(t) 3 * t,
                              function(t) log1p(t),
                              function(t) -4 * t + 1,
                              function(t) 5 * t,
                              function(t) 2 * t^2,
                              function(t) sin(pi * t))){
  # Time points are a sequence between 0 and 1
  time <- seq(0, maxt, 0.05)

  # Simulation parameters for setting Y0
  ranges <- range_start:range_end

  # How many missing times?
  miss_n <- seq(miss_p, 1, by = 0.1)

  # Set up the longitudinal correlation structure
  # What are these numbers?
  cor_matrix <- cor_user(length(time),
                         user_var = user_var,
                         cor_str = cor_str,
                         rho = rho)

  # Simulate data.

  sim_data_missing <- list()

  # Generate data for all i.
  for (i in 1:n) {
    sim_data <- NULL


    # Simulate external variable Z
    if (Z_type == "binary") {
      Z <- stats::rbinom(n = length(time),
                         size = 1,
                         prob = prob1)
    } else if (Z_type == "continuous") {
      Z <- stats::runif(n = length(time), min = 0, max = 3)
    } else {
      stop("Z_type must be either binary or continuous")
    }

    # Generate alphas based on the baseline and slope functions specified
    alpha <- data.frame(matrix(ncol = K, nrow = length(time)))
    for (k in 1:K) {
      baseline_value <- baseline_fxns[[k]](time)
      slope_value  <- slope_fxns[[k]](time)
      alpha[,k] <- exp(baseline_value + Z * slope_value)
    }
    names(alpha) <- paste0("a", 1:K)

    # Get the error matrix representing the longitudinal correlation
    longitudinal_cor_error <- t(MASS::mvrnorm(n = K,
                                              mu = rep(0, length(time)),
                                              Sigma = cor_matrix,
                                              tol = 1e-6,
                                              empirical = FALSE,
                                              EISPACK = FALSE))

    # Generate error structure for each j: NOTE: this means that there can't be
    # simulated error structures where j != j'
    for (j in 1:length(time)) {
      if (length(ranges) == 1) {
        Y0 <- ranges
      } else {
        Y0 <- sample(ranges, 1)
      }
      # Generate the compositional correlation
      dm_counts <- generate_DM_counts(Y0, as.numeric(alpha[j, ]))

      # Add the longitudinal correlation to the compositional correlation
      # Round it to be a count
      dm_witherror <- round(dm_counts + longitudinal_cor_error[j, ], 0)

      # Make sure no zeros.
      dm_witherror[dm_witherror < 0] <- 0

      # Define the total sum
      #total_n <- rowSums(dm_witherror)

      # Return data with the total n and id.
      dm_new <- c(i, Z[j], dm_witherror)
      sim_data <- rbind(sim_data, dm_new)
    }

    # Include missingness.
    # Set missingness
    miss_rate <- sample(miss_n, 1)
    newtime <- sort(sample(1:length(time), ceiling(miss_rate * length(time))))
    sim_data_missing[[i]] <- data.frame(time, sim_data)[newtime, ]

    rownames(sim_data_missing[[i]]) <- NULL
  }

  # Organize data.
  data <- dplyr::bind_rows(sim_data_missing)
  names(data) <- c("time",
                   "individual",
                   "Z",
                   paste0("Response_", 1:K))

  #data$Capture.Number <- data$time * (length(unique(data$time)) - 1) + 1
  data$individual <- as.factor(data$individual)

  Y <- dplyr::select(data, -c(.data$Z))
  Z <- data$Z
  return(list(Y = Y,
              Z = Z))
}







#  Simulation helpers: ---------
#' Simulate Dirichlet multinomial counts.
#'
#' @param Y0 Total sum of counts in DM
#' @param alpha Vector of Dirichlet parameters
#'
#' @returns DM counts
#' @export
generate_DM_counts <- function(Y0, alpha) {
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
#' @param mi number of time points for ith sample
#' @param user_var User supplied variance
#' @param rho rho correlation parameter
#' @param cor_str specified correaltion structure (of IND, CON-d, or AR1d)
#'
#' @returns User defined correlation matrix (for a single sample i)
#' @export
cor_user <- function(mi, user_var, cor_str, rho) {
  if (cor_str == "CON-d") {
    identity_matrix <- diag(mi)
    ones_matrix <- matrix(1, nrow = mi, ncol = mi)
    cor_structure <- (ones_matrix - identity_matrix) * rho + identity_matrix
  } else if (cor_str == "AR1-d") {
    row_indices <- matrix(rep(1:mi, each = mi), nrow = mi)
    col_indices <- matrix(rep(1:mi, each = mi), ncol = mi, byrow = TRUE)
    # Compute the absolute difference between indeces: |j1-j2|
    distance_matrix <- abs(row_indices - col_indices)
    # returns (rho^|j1-j2|)
    cor_structure <- rho^distance_matrix
  }
  else if (cor_str == "IND") {
    identity_matrix <- diag(mi)
    cor_structure <- identity_matrix
  } else
    stop("Invalid cor_str")

  # Scale by user_var
  cor <- user_var * cor_structure
  return(cor)
}



# Diagnostics  ------------------


#' Look at clustering results if true clusters are known
#'
#' @param res model object from cleverly
#' @param true_cluster vector with the true cluster membership
#'
#' @returns res model with additional list item diagnostics that contains clustering results of rand, adj.rand, jaccard, CER, and nclust
#' @export
get_cluster_diagnostics <- function(res, true_cluster){

  found_cluster <- res$cluster$membership

  rand <- fossil::rand.index(found_cluster, true_cluster)
  adj_rand <- mclust::adjustedRandIndex(found_cluster, true_cluster)
  jaccard <- length(intersect(found_cluster, true_cluster)) / length(union(found_cluster, true_cluster))
  CER <- mclust::classError(classification = found_cluster,
                            class = true_cluster)$errorRate

  res$cluster_diagnostics <- data.frame("rand" = rand,
                                        "adj_rand" = adj_rand,
                                        "jaccard" = jaccard,
                                        "CER" = CER,
                                        "nclust" = res$cluster$no)
  return(res)
}


#' Get proportion correct cluster possible.
#'
#' @param res cleverly object
#' @param true_cluster true cluster membership
#'
#' @returns Proportion of simulation runs the correct cluster was an option
#' @export
correct_cluster_option <- function(res, true_cluster = rep(1:3, each = 4)){
  val <- purrr::map(res, "possible_clusters") %>%
    purrr::map_dfc(~map(.x, ~.x$membership == true_cluster ) %>%
              map_dbl(~sum(.x)==12)) %>%
    colSums() %>% unname() %>% suppressMessages()

  return(mean(val > 0))
}


#' get performance summary table
#'
#' @param list simulation list
#'
#' @returns Returns the object with the performance summary table
#' @export
get_performance_summary <- function(list){

  table <- purrr::map_dfr(list, ~map_dfr(.x,"cluster_diagnostics") %>%
                            summarize_all(list("mean" = mean, "sd" = sd)),
                          .id = "corstr")


  return(table)

}


# Visualizations ----------------------------------------------------------

#' Plot simulation data
#'
#' @param sim simulated data
#' @param K number of responses
#'
#' @returns ggplot object with the simulated data for each response
#' @export

plot_sim_data <- function(sim, K = 12){

  binary_Z <- length(unique(sim$Z)) <= 2

  if (binary_Z) {
    plot <- sim %>%
      tidyr::pivot_longer(-c(.data$individual,
                             .data$time,
                             .data$Capture.Number,
                             .data$total_n,
                             .data$Z)) %>%
      dplyr::mutate(name = factor(.data$name,
                                  labels = paste0("Response ", 1:K))) %>%
      ggplot2::ggplot(ggplot2::aes(x = time,
                                   y = value,
                                   color = factor(Z))) +
      ggplot2::geom_jitter(size = 1) +
      ggplot2::facet_wrap(~name) +
      ggplot2::labs(title = "Simulated Data",
                    color = "Z",
                    shape = "Z",
                    y = "Count",
                    x = "Time")
  } else {
    plot <- sim %>%
      tidyr::pivot_longer(-c(.data$individual,
                             .data$time,
                             .data$Capture.Number,
                             .data$total_n,
                             .data$Z)) %>%
      dplyr::mutate(name = factor(.data$name,
                                  labels = paste0("Response ", 1:K))) %>%
      ggplot2::ggplot(ggplot2::aes(x = time,
                                   y = value,
                                   color = Z)) +
      ggplot2::geom_jitter(size = 1) +
      ggplot2::facet_wrap(~name) +
      ggplot2::labs(title = "Simulated Data",
                    color = "Z",
                    shape = "Z",
                    y = "Count",
                    x = "Time")
  }

  return(plot)
}
