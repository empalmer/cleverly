
test_that("Gradient il dimension",{
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  P <- 6
  mi_vec <- rep(3, 5)
  i <- 1
  i_index <- 1

  Y0 <- rowSums(Y)
  alpha <- get_alpha_list(beta, Z, B, K, i_index = 1, mi_vec = mi_vec, L = 2, P = P)

  partials_il <- get_partials_il(i = i,
                                 l = 0,
                                 Y0,
                                 Z,
                                 B,
                                 beta,
                                 alpha_i = alpha[[i]],
                                 mi_vec,
                                 i_index,
                                 P,
                                 K)

  gradient_il <- get_gradient_il(i = 1,
                                 l = 2,
                                 Y = Y,
                                 Y0 = Y0,
                                 mi_vec = mi_vec,
                                 i_index = c(0, cumsum(mi_vec)),
                                 phi = .5,
                                 beta = beta,
                                 alpha = alpha,
                                 Z = Z,
                                 B = B,
                                 Vi_inv = diag(1, 12),
                                 partials_il = partials_il,
                                 K = K)

  # Check dimensions
  expect_length(gradient_il, P*K)
})


