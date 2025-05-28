
#' CLR cluster
#'
#' Note: only implemented for the case of 1 external variable.
#'
#' @param y
#' @param Z
#' @param beta
#' @param lp
#' @param lp_minus
#' @param B
#' @param clusters
#' @param K
#' @param P
#' @param M
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
  # Build design matrix for each outcome: [A | B], size M x 12
  X <- kronecker(diag(K), ZB)  # (M*K) x (12*p)

  # Stack responses
  Y_vec <- as.vector(y_clr)

  C <- matrix(0, nrow = (L + 1) * P, ncol = (L + 1) * P * K)

  # Fill in C for beta_0 constraints (rows 1:6)
  for (j in 1:P) {
    for (k in 1:K) {
      idx <- (k - 1) * (L + 1) * P + j  # location of beta_0[j, k]
      C[j, idx] <- 1
    }
  }
  # Fill in C for beta_1 constraints (rows 7:12)
  for (j in 1:P) {
    for (k in 1:K) {
      idx <- (k - 1) * (L + 1) * P + P + j  # location of beta_1[j, k]
      C[P + j, idx] <- 1
    }
  }

  # KKT system
  XtX <- t(X) %*% X
  XtY <- t(X) %*% Y_vec
  zero_block <- matrix(0, nrow = nrow(C), ncol = nrow(C))

  KKT_mat <- rbind(
    cbind(XtX, t(C)),
    cbind(C, zero_block)
  )
  rhs <- rbind(XtY, matrix(0, nrow = nrow(C), ncol = 1))
  sol <- solve(KKT_mat, rhs)
  beta_constrained <- sol[1:((L + 1) * P * K), ]


  # Reshape into 2 coefficient matrices
  beta0 <- matrix(NA, nrow = P, ncol = K)
  beta1 <- matrix(NA, nrow = P, ncol = K)
  for (k in 1:K) {
    beta_k <- beta_constrained[((k - 1) * (L + 1) * P + 1):(k * (L + 1) * P)]
    beta0[, k] <- beta_k[1:P]
    beta1[, k] <- beta_k[(P + 1):((L + 1) * P)]
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
    dplyr::mutate(response = factor(response, levels = 1:K))

  # Compute fitted baseline values.
  y_hat_baseline <- matrix(NA, nrow = nrow(Y), ncol = K)
  for (k in 1:K) {
    y_hat_baseline[, k] <- ZB_list[[1]] %*% beta0[, k]
  }
  y_hat_baseline_df <- data.frame(time = time, Z = Z[, 2], y_hat_baseline, y_clr)
  colnames(y_hat_baseline_df) <- c("time",
                                   "Z",
                                   paste0("yhat_", 1:K),
                                   paste0("y_", 1:K))
  y_hat_baseline_df <- y_hat_baseline_df %>%
    tidyr::pivot_longer(
      cols = tidyr::matches("yhat|y"),  # Selects both yhat and y columns
      names_to = c(".value", "response"),
      names_pattern = "(yhat|y)_(\\d+)"  # Splits into two parts: yhat/y and the number
    ) %>%
    dplyr::mutate(response = factor(response, levels = 1:K))





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



#' CLR cluster
#'
#' @param y
#' @param Z
#' @param beta
#' @param lp
#' @param lp_minus
#' @param B
#' @param clusters
#' @param K
#' @param P
#' @param M
#'
#' @returns
#' @export
logY_ols_cluster <- function(Y, Z, lp,  lp_minus, B, K, P, M) {

  L <- ncol(Z) - 1


  gamma <- matrix(nrow = P * (L + 1), ncol = K)


  B <- get_B(time = sim$time,
             order = 3,
             nknots = 3)

  ZB_list <- list()
  for (l in 0:L) {
    Z_l <- diag(Z[,l + 1])
    ZB_list[[l + 1]] <- Z_l %*% B
  }
  ZB <- do.call(cbind, ZB_list)

  Y <- dplyr::select(Y,-c(time, individual))
  logY <- log(Y + 1)

  for (k in 1:K) {

  }

}
