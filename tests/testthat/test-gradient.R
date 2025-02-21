test_that("Yi  minus mui", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)

  term <- get_Yi_minus_mui(i = 1,
                           Y = Y,
                           mi_vec = mi_vec,
                           beta = beta,
                           Z = Z,
                           B = B,
                           K = K)
  expect_length(term, K*mi)
})


test_that("Yi  minus mui SIM", {
  sim <- base_sim()

  i <- 3
  term <- get_Yi_minus_mui(i = i,
                           Y = sim$Y,
                           mi_vec = sim$mi_vec,
                           beta = sim$beta,
                           Z = sim$Z,
                           B = sim$B,
                           K = sim$K)

  expect_length(term, sim$K*sim$mi_vec[i])
})



test_that("Vi inverse", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)

  V_i_inv <- get_Vi_inv(i = 1,
                        Y = Y,
                        mi_vec = mi_vec,
                        phi = .5,
                        beta = beta,
                        Z = Z,
                        B = B,
                        K = K)

  # Check dimensions
  expect_equal(dim(V_i_inv)[1], K*mi)
})


test_that("Vi inv SIM ", {
  sim <- base_sim()

  i <- 1
  V_i_inv <- get_Vi_inv(i = i,
                        Y = sim$Y,
                        mi_vec = sim$mi_vec,
                        phi = .5,
                        beta = sim$beta,
                        Z = sim$Z,
                        B = sim$B,
                        K = sim$K)

  # Check dimensions
  expect_equal(dim(V_i_inv)[1], sim$K*sim$mi_vec[i])
})



test_that("Partials ijl dimension", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6
  mi_vec <- rep(3, 5)

  partials <- get_partials_ijl(i = 1,
                           j = 2,
                           l = 2,
                           Y = Y,
                           mi_vec = mi_vec,
                           beta = beta,
                           Z = Z,
                           B = B,
                           K = K)

  # Check dimensions
  expect_equal(dim(partials)[1], K*P)
  expect_equal(dim(partials)[2], K)
})

test_that("Partials ijl dim SIM", {
  sim <- base_sim()
  partials <- get_partials_ijl(i = 1,
                               j = 2,
                               l = 2,
                               Y = sim$Y,
                               mi_vec = sim$mi_vec,
                               beta = sim$beta,
                               Z = sim$Z,
                               B = sim$B,
                               K = sim$K)

  # Check dimensions
  expect_equal(dim(partials)[1], sim$K*sim$P)
  expect_equal(dim(partials)[2], sim$K)
})



test_that("Partials il dimension",{
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6
  mi_vec <- rep(3, 5)

  partials <- get_partials_il(i = 1,
                              l = 2,
                              Y = Y,
                              mi_vec = mi_vec,
                              beta = beta,
                              Z = Z,
                              B = B)

  # Check dimensions
  expect_equal(dim(partials)[1], K*P)
  expect_equal(dim(partials)[2], K*mi)
})


test_that("Partials il dim SIM ",{
  sim <- base_sim()

  i <- 1
  partials <- get_partials_il(i = i,
                              l = 2,
                              Y = sim$Y,
                              mi_vec = sim$mi_vec,
                              beta = sim$beta,
                              Z = sim$Z,
                              B = sim$B)

  # Check dimensions
  expect_equal(dim(partials)[1], sim$K*sim$P)
  expect_equal(dim(partials)[2], sim$K*sim$mi_vec[i])
})





test_that("Gradient il dimension",{
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6
  mi_vec <- rep(3, 5)

  gradient_il <- get_gradient_il(i = 1,
                              l = 2,
                              Y = Y,
                              mi_vec = mi_vec,
                              phi = .5,
                              beta = beta,
                              Z = Z,
                              B = B)

  # Check dimensions
  expect_length(gradient_il, P*K)
})


test_that("Gradient il dim SIM",{
  sim <- base_sim()

  gradient_il <- get_gradient_il(i = 1,
                                 l = 2,
                                 Y = sim$Y,
                                 mi_vec = sim$mi_vec,
                                 phi = .5,
                                 beta = sim$beta,
                                 Z = sim$Z,
                                 B = sim$B)

  # Check dimensions
  expect_length(gradient_il, sim$P*sim$K)
})



test_that("Gradient i",{
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6
  mi_vec <- rep(3, 5)

  gradient_il <- get_gradient_il(i = 1,
                                 l = 2,
                                 Y = Y,
                                 mi_vec = mi_vec,
                                 phi = .5,
                                 beta = beta,
                                 Z = Z,
                                 B = B)

  # Check dimensions
  expect_length(gradient_il, P*K)
})


test_that("Gradient l",{
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6
  sid <- rep(1:5, each = 3)
  tid <- rep(1:3, 5)


  gradient_l <- get_gradient_l(Y = Y,
                               mi_vec = rep(3, 5),
                               l = 0,
                               phi = .5,
                               beta = beta,
                               Z = Z,
                               B = B)

  # Check dimensions
  expect_length(gradient_l, P*K)
})

test_that("Gradient l SIM",{
  sim <- base_sim()

  gradient_l <- get_gradient_l(Y = sim$Y,
                               mi_vec = rep(3, 5),
                               l = 0,
                               phi = .5,
                               beta = sim$beta,
                               Z = sim$Z,
                               B = sim$B)

  # Check dimensions
  expect_length(gradient_l, sim$P*sim$K)
})
