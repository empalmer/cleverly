test_that("Yi  minus mui", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3

  term <- get_Yi_minus_mui(i = 1,
                           Y = Y,
                           mi = mi,
                           beta = beta,
                           Z = Z,
                           B = B,
                           K = K)
  expect_length(term, K*mi)
})


test_that("Vi inverse", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3


  V_i_inv <- get_Vi_inv(i = 1,
                        Y = Y,
                        mi = mi,
                        phi = .5,
                        beta = beta,
                        Z = Z,
                        B = B,
                        K = K)

  # Check dimensions
  expect_equal(dim(V_i_inv)[1], K*mi)
})



test_that("Partials ijl dimension", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6

  partials <- get_partials_ijl(i = 1,
                           j = 2,
                           l = 2,
                           Y = Y,
                           mi = mi,
                           beta = beta,
                           Z = Z,
                           B = B,
                           K = K)

  # Check dimensions
  expect_equal(dim(partials)[1], K*P)
  expect_equal(dim(partials)[2], K)
})


test_that("Partials il dimension",{
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6

  partials <- get_partials_il(i = 1,
                              l = 2,
                              Y = Y,
                              mi = mi,
                              beta = beta,
                              Z = Z,
                              B = B)

  # Check dimensions
  expect_equal(dim(partials)[1], K*P)
  expect_equal(dim(partials)[2], K*mi)
})

test_that("Gradient il dimension",{
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6

  gradient_il <- get_gradient_il(i = 1,
                              l = 2,
                              Y = Y,
                              mi = mi,
                              phi = .5,
                              beta = beta,
                              Z = Z,
                              B = B)

  # Check dimensions
  expect_length(gradient_il, P*K)
})

test_that("Gradient i",{
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6

  gradient_il <- get_gradient_il(i = 1,
                                 l = 2,
                                 Y = Y,
                                 mi = mi,
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
                               subject_ids = sid,
                               time_ids = tid,
                               l = 0,
                               phi = .5,
                               beta = beta,
                               Z = Z,
                               B = B)

  # Check dimensions
  expect_length(gradient_l, P*K)
})

