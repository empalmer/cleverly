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
