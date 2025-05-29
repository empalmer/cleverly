
#' CLR cluster
#'
#' Note: only implemented for the case of 1 external variable.
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param Z Full Z matrix (including intercept)
#' @param lp Which column is clustering can only be 0 or 1
#' @param B B spline basis matrix
#' @param K Number of responses
#' @param P number of B-spline coefficients (order + nknots)
#' @param M Number of samples times time points for each sample
#' @param time Vector of numeric time values.
#'
#' @returns
#' @export
CLR_cluster <- function(Y, Z, time, B, lp, K, P, M) {

  # Get design matrix ready:
  L <- ncol(Z) - 1
  ZB_list <- list()
  for (l in 0:L) {
    Z_l <- diag(Z[,l + 1])
    ZB_list[[l + 1]] <- Z_l %*% B
  }
  ZB <- do.call(cbind, ZB_list)

  # Transform Y:
  y_clr <- compositions::clr(Y)

  # Multivariate constrained optimization
  X <- kronecker(diag(K), ZB)  # (M*K) x (12*p)

  # Vectorize y (not multivariate)
  Y_vec <- as.vector(y_clr)

  # Make constraint matrix
  E <- kronecker(matrix(1, ncol = K), diag(K))


  # Lsei solves min(|AX - B|2) subject to EX = F
  lsei <- limSolve::lsei(A = X,
                         B = Y_vec,
                         E = E,
                         F = rep(0, nrow(E)))


  # lsei outputs one long beta vector, so we need to reshape it
  beta_constrained <- lsei$X
  beta0 <- matrix(NA, nrow = P, ncol = K)
  beta1 <- matrix(NA, nrow = P, ncol = K)

  for (k in 1:K) {
    # Indices for beta_0 and beta_1 for outcome k
    idx_start <- (k - 1) * 2 * P + 1
    idx_end <- k * 2 * P
    beta_k <- beta_constrained[idx_start:idx_end]

    beta0[, k] <- beta_k[1:P]             # first P values = beta_0
    beta1[, k] <- beta_k[(P + 1):(2 * P)] # next P values = beta_1
  }


  # Compute fitted values
  y_hat <- matrix(NA, nrow = nrow(Y), ncol = K)
  for (k in 1:K) {
    y_hat[, k] <- ZB_list[[1]] %*% beta0[, k] + ZB_list[[2]] %*% beta1[, k]
  }
  y_hat_df <- data.frame(time = time, Z = Z[, 2], y_hat, y_clr)
  colnames(y_hat_df) <- c("time",
                          "Z",
                          paste0("yhat_", 1:K),
                          paste0("y_", 1:K))
  y_hat_df <- y_hat_df %>%
    tidyr::pivot_longer(
      cols = tidyr::matches("yhat|y"),  # Selects both yhat and y columns
      names_to = c(".value", "response"),
      names_pattern = "(yhat|y)_(\\d+)"  # Splits into two parts: yhat/y and the number
    ) %>%
    dplyr::mutate(response = factor(.data$response, levels = 1:K))

  # Compute fitted baseline values.
  y_hat_clr_baseline <- matrix(NA, nrow = nrow(Y), ncol = K)
  for (k in 1:K) {
    y_hat_clr_baseline[, k] <- ZB_list[[1]] %*% beta0[, k]
  }

  # Transform back into relative abundances.
  y_hat_baseline <- compositions::clrInv(y_hat_clr_baseline)
  #y <- compositions::clrInv(y_clr)
  y <- Y/rowSums(Y)

  y_hat_baseline_df <- data.frame(time = time, Z = Z[, 2], y_hat_baseline, y, y_hat_clr_baseline, y_clr)

  colnames(y_hat_baseline_df) <- c("time",
                                   "Z",
                                   paste0("yhat_", 1:K),
                                   paste0("y_", 1:K),
                                   paste0("yhatclr_", 1:K),
                                   paste0("yclr_", 1:K))
  y_hat_baseline_df <- y_hat_baseline_df %>%
    tidyr::pivot_longer(
      cols = tidyr::matches("yhat|y"),  # Selects both yhat and y columns
      names_to = c(".value", "response"),
      names_pattern = "(yhat|y|yhatclr|yclr)_(\\d+)"  # Splits into two parts: yhat/y and the number
    ) %>%
    dplyr::mutate(response = factor(.data$response, levels = 1:K))

  # Perform k means clustering using gap statistic to determine number of clusters
  if (lp == 0) {
    beta <- t(beta0)
  } else {
    beta <- t(beta1)
  }

  gap_stat <- cluster::clusGap(beta, FUN = kmeans, nstart = 25, K.max = K - 1)
  no <- cluster::maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
  kmeans_res <- kmeans(beta, centers = no, nstart = 25)


  res <- list(clusters = list(membership = kmeans_res$cluster,
                              no = no),
              y = y_clr,
              y_hat = y_hat_df,
              y_hat_baseline = y_hat_baseline_df,
              beta0 = beta0,
              beta1 = beta1)

  return(res)
}

