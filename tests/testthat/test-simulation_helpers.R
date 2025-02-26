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


test_that("No Z - Chenyangs code", {

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





  res <- cleverly(Y = Y,
                  subject_ids = id,
                  time = time,
                  gammas = 1,
                  psi = 800,
                  phi = 1,
                  max_admm_iter = 5,
                  max_outer_iter = 1)

  est_beta <- res$result$beta
  phis <- res$result$phis_list[[1]]
  v <- res$result$v

  # View the betas to see if they are clustered.
  data.frame(est_beta) %>%
    dplyr::mutate(K = rep(LETTERS[1:12],each = 6),
                  id = rep(1:6, 12)) %>%
    tidyr::pivot_wider(names_from = K,
                       values_from = est_beta)

  y_hat <- res$result$y_hat
  colnames(y_hat) <- c("Taxa.1","Taxa.2",
                       "Taxa.3","Taxa.4",
                       "Taxa.5","Taxa.6",
                       "Taxa.7","Taxa.8",
                       "Taxa.9","Taxa.10",
                       "Taxa.11","Taxa.12")

  #Un-normalize?
  y_ra <- Y/rowSums(Y)
  colnames(y_ra) <- c("Taxa.1","Taxa.2",
                       "Taxa.3","Taxa.4",
                       "Taxa.5","Taxa.6",
                       "Taxa.7","Taxa.8",
                       "Taxa.9","Taxa.10",
                       "Taxa.11","Taxa.12")
  y_ra <- data.frame(y_ra) %>%
    dplyr::mutate(time = time) %>%
    tidyr::pivot_longer(-time)



  y_hat_ra <- data.frame(y_hat) %>%
    dplyr::mutate(time = time) %>%
    tidyr::pivot_longer(-time)


  # Y-hats
  Y_ra_df <- data.frame(y_ra,
                        y_hat = y_hat_ra$value) %>%
    dplyr::mutate(name = factor(name,
                                levels = c("Taxa.1","Taxa.2",
                                           "Taxa.3","Taxa.4",
                                           "Taxa.5","Taxa.6",
                                           "Taxa.7","Taxa.8",
                                           "Taxa.9","Taxa.10",
                                           "Taxa.11","Taxa.12")))
  Y_ra_df %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = time, y = value),
                        size = 1, color = "black") +
    ggplot2::geom_line(ggplot2::aes(x = time, y = y_hat),
                       linewidth = 2, color = "blue") +
    ggplot2::facet_wrap(~name) +
    ggplot2::labs(title = "After 2 iterations",
                  x = "Time",
                  y = "ra")


})
