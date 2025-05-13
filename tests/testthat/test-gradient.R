test_that("Yi  minus mui", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  Y0 <- colSums(Y)
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  L <- 1
  P <- 6
  alpha <- get_alpha_list(beta,
                          Z,
                          B,
                          K,
                          i_index,
                          mi_vec,
                          L, P)

  term <- get_Yi_minus_mui(i = 1,
                           Y = Y,
                           Y0 = Y0,
                           mi_vec = mi_vec,
                           i_index = i_index,
                           alpha = alpha,
                           K = K)
  expect_length(term, K * mi)
})




test_that("Vi inverse", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  V_i_inv <- get_Vi_inv(i = 1,
                        Y = Y,
                        mi_vec = mi_vec,
                        i_index = i_index,
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
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  partials <- get_partials_ijl(i = 1,
                           j = 2,
                           l = 2,
                           Y = Y,
                           mi_vec = mi_vec,
                           i_index = i_index,
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
  mi_vec <- rep(3, 5)

  partials <- get_partials_il(i = 1,
                              l = 2,
                              Y = Y,
                              mi_vec = mi_vec,
                              i_index = c(0, cumsum(mi_vec)),
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
  mi_vec <- rep(3, 5)

  gradient_il <- get_gradient_il(i = 1,
                              l = 2,
                              Y = Y,
                              mi_vec = mi_vec,
                              i_index = c(0, cumsum(mi_vec)),
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
  Y0 <- colSums(Y)
  K <- 4
  mi <- 3
  P <- 6
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  gradient_il <- get_gradient_il(i = 1,
                                 l = 2,
                                 Y = Y,
                                 Y0 = Y0,
                                 K = K,
                                 mi_vec = mi_vec,
                                 i_index = i_index,
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
  Y0 <- colSums(Y)
  K <- 4
  mi <- 3
  P <- 6
  sid <- rep(1:5, each = 3)
  tid <- rep(1:3, 5)
  i_index <- 1
  mi_vec <- rep(3, 5)
  L <- 1
  alpha <- get_alpha_list(beta, Z, B, K, i_index, mi_vec, L, P)
  l = 0
  partials_l <- get_partials_l_list(Y0,
                                    l,
                                    mi_vec,
                                    i_index,
                                    beta,
                                    alpha,
                                    Z,
                                    B,
                                    P,
                                    K)

  gradient_l <- get_gradient_l(Y = Y,
                               Y0 = Y0,
                               mi_vec = rep(3, 5),
                               i_index = c(0, cumsum(rep(3, 5))),
                               l = 0,
                               alpha = alpha,
                               partials_l = partials_l,
                               phi = .5,
                               beta = beta,
                               K = K,
                               Z = Z,
                               B = B,
                               P = P)

  # Check dimensions
  expect_length(gradient_l, P*K)
})
