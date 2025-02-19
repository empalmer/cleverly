test_that("timepoints work", {
  time <- sim_timepoints(n = 5)

  expect_type(time, "list")
})


test_that("Z works", {
  time <- sim_timepoints(n = 5)

  mis <- time$mis

  Z <- sim_Z(mis)
})


test_that("sim_Y_ij works", {
  time_list <- sim_timepoints(n = 5)
  mis <- time_list$mis
  time <- time_list$X$time

  K <- 4
  Z <- sim_Z(mis)

  B <- get_B(time, order = 3, nknots = 3)

  # 2 clusters
  betaC1 <- matrix(c(rep(c(1, 1, 1, 1, 1, 1), 2),     #l = 0
                     rep(c(-2, -2, -2, -2, -2, -2),2),  #l = 1
                     rep(c(1, 2, 3, 4, 5, 6), 2)), ncol = 3)  #l = 2
  betaC2 <- -betaC1

  beta <- rbind(betaC1, betaC2)

  Y_ij <- sim_Y_ij(i = 1,
                   j = 1,
                   beta,
                   Z,
                   B,
                   Y_ij0 = 100,
                   K = K,
                   mi = mis[1])

  Y_i <- sim_Yi(i = 1,
                beta = beta,
                Z = Z,
                B = B,
                K = K,
                mi = mis[1])

  Y <- sim_Y(beta = beta,
             Z = Z,
             B = B,
             K = K,
             mis = mis)
})
