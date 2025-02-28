test_that("timepoints work", {
  time <- sim_timepoints(n = 5)

  expect_type(time, "list")
})


test_that("Z works", {
  time <- sim_timepoints(n = 5)

  mi_vec <- time$mi_vec

  Z_sim <- sim_Z(mi_vec)
})


test_that("Test sim Y works", {
  set.seed(124)
  time_list <- sim_timepoints(n = 5)
  mi_vec <- time_list$mi_vec
  time <- time_list$X$time

  K <- 4
  Z <- sim_Z(mi_vec)

  B <- get_B(time, order = 3, nknots = 3)

  # 2 clusters
  betaC1 <- matrix(c(rep(c(1, 1, 1, 1, 1, 1), 2),     #l = 0
                     rep(c(-2, -2, -2, -2, -2, -2),2),  #l = 1
                     rep(c(1, 2, 3, 4, 5, 6), 2)), ncol = 3)  #l = 2
  betaC2 <- .5*betaC1

  beta <- rbind(betaC1, betaC2)

  Y_ij <- sim_Y_ij(i = 1,
                   j = 1,
                   beta,
                   Z,
                   B,
                   Y_ij0 = 100,
                   K = K,
                   mi_vec = mi_vec)

  Y_i <- sim_Yi(i = 1,
                beta = beta,
                Z = Z,
                B = B,
                K = K,
                mi_vec = mi_vec)

  Y <- sim_Y(beta = beta,
             Z = Z,
             B = B,
             K = K,
             mi_vec = mi_vec)
})


test_that("No Z - Functional simulation", {
  testthat::skip("Debugging, not test")

  # Generate simulation data
  set.seed(123)
  sim <- sim_noZ()
  time <- sim$time
  id <- sim$individual
  time_id <- sim$Capture.Number
  Yij0 <- sim$total_n
  Y <- dplyr::select(sim, -c(
                        "individual",
                        "time",
                        "total_n",
                        "Capture.Number"))
  Y <- as.matrix(Y)
  K <- ncol(Y)

  # Run algorithm
  iter <- 25
  res <- cleverly(Y = Y,
                  subject_ids = id,
                  time = time,
                  gammas = 1000, # controls smoothness
                  tau = 1/2,
                  theta = 300,
                  psi = 8000, # controls clustering
                  C = 100,
                  max_admm_iter = iter,
                  max_outer_iter = 1)


  # Diagnostics:
  phis <- res$result$phis_list[[1]]
  v <- res$result$v
  v
  v_mat <- matrix(v, nrow = 6)

  # See the final fit
  visualize_final_fit(Y = Y,
                      res = res)

  # View how curves change over iterations.
  betas <- res$result$admm_beta_list[[1]]
  beta_path(betas,
            K = K,
            time = time)

  # Examine differences:
  # Differences in the ADMM iterations
  diffs <- res$result$admm_diffs[[1]]
  plot_differneces(diffs)
  diffs

  # Examine clusters:
  A <- res$result$clusters

})



test_that("Bspline sim no Z",{
  testthat::skip("Debugging, not test")
  set.seed(123)
  n <- 100
  time_list <- sim_timepoints(n = n)
  time <- time_list$X$time
  mi_vec <- time_list$mi_vec
  is <- rep(1:n, mi_vec)
  M <- sum(mi_vec)
  K <- 4
  Z <- matrix(rep(1, M), nrow = M)
  B <- get_B(time, order = 3, nknots = 3)
  P <- 6

  # 2 clusters
  betaC1 <- matrix(c(rep(c(1, 2, 3, 4, 5, 6), 2), ncol = 1))#l = 2
  betaC2 <- matrix(.5*rev(betaC1), ncol = 1)
  beta <- rbind(betaC1, betaC2)

  Y <- sim_Y(beta = beta,
             Z = Z,
             B = B,
             K = K,
             mi_vec = mi_vec)
  Y_ra <- Y/rowSums(Y)

  # Plot data:
  data.frame(time = time, Y_ra) %>%
    tidyr::pivot_longer(-time) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~name)

  # Truth:
  y_hat_truth <- estimate_y(beta = beta,
                       B = B,
                       Z = Z,
                       K = K)
  data.frame(time = time, y_hat_truth) %>%
    tidyr::pivot_longer(-time) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~name)


  # Run code
  iter <- 25
  res <- cleverly(Y = Y,
                  subject_ids = is,
                  time = time,
                  gammas = .1,
                  psi = .01,
                  max_admm_iter = iter,
                  max_outer_iter = 1)

  visualize_curve(beta_new, B, Z, K = 4, time)

  betas <- res$result$admm_beta_list[[1]]
  plot <- beta_path(betas, K)


})
