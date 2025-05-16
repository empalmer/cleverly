test_that("Algorithm 1", {
  skip("TODO - alg1")
  sim <- simulation_data(n = 2,
                         range_start = 5000,
                         range_end = 6000,
                         nknots = 3,
                         K = 2,
                         order = 3,
                         maxt = .15,
                         user_var = 100,
                         cor_str = "CON-d",
                         rho = 0.5,
                         prob1 = .5,
                         baseline_fxns = list(
                           function(t) t,
                           function(t) -3 * t + 3
                         ),
                         slope_fxns = list(
                           function(t) 2,
                           function(t) t
                         ))
  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))
  Z <- format_Z(sim$Z)

  individual <- as.numeric(sim$individual)
  time <- sim$time
  mi_vec <- get_mi_vec(Y, individual, time)$mi
  result <- algorithm1(Y = Y,
                        Z = Z,
                        time = time,
                        mi_vec = mi_vec,
                        lp = 0,
                        cor_str = "IND",
                        gammas = c(1,1),
                        psi = 800,
                        tau = 0.005,
                        theta = 300,
                        C = 100,
                        d = 2,
                        nknots = 3,
                        order = 3,
                        epsilon_b = 1,
                        epsilon_r = 1,
                        epsilon_d = 1,
                        epsilon_2 = 1,
                        run_min = 1,
                        max_outer_iter = 1,
                        max_admm_iter = 1,
                        max_2_iter = 1)

  expect_true(is.list(result))

})


# Initializing algorithm  -----------------------------------------------


test_that("format Z", {
  Z_w1 <- readRDS(test_path("test_data", "Z.rds"))

  # Case when Z includes column of 1s
  Z <- format_Z(Z_w1)
  expect_equal(Z, Z_w1)

  # Case when Z does not include a column of 1s
  Z_wo1 <- Z_w1[,-1]
  Z1 <- format_Z(Z_wo1)
  expect_equal(Z1, Z_w1)

})


test_that("D matrix", {
  K <- 5
  D <- get_D(K = K, d = 2, order = 3, nknots = 3)

  P <- 6
  expect_equal(dim(D), c(P * K, P*K))
})


test_that("A matrix", {
  K <- 4
  Kappa <- t(combn(K,2))

  P <- 6
  Kappa_size <- nrow(Kappa)
  A <- get_A(Kappa, K, P)

  expect_equal(dim(A), c(nrow(Kappa)*P, P*K))
})
