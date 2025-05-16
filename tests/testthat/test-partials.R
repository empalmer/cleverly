test_that("get_partials_ijl returns matrix with correct dimensions", {
  i <- 1; j <- 1; l <- 0
  mi_vec <- c(2)
  i_index <- c(0)
  Y0 <- c(100, 150)
  Z <- matrix(1, nrow = 2, ncol = 2)
  B <- matrix(1:6, nrow = 2, ncol = 3)
  P <- 3; K <- 2
  beta <- matrix(1:(P*K), nrow = P*K, ncol = 1)
  alpha_ij <- rep(1, K)

  partials <- get_partials_ijl(i, j, l, mi_vec, i_index, Y0, Z, B, beta, alpha_ij)
  expect_type(partials, "double")
  expect_equal(dim(partials), c(P*K, K))
})


test_that("get_partials_il returns matrix with correct dimensions", {
  i <- 1; l <- 0
  mi_vec <- c(2)
  i_index <- c(0)
  Y0 <- c(100, 150)
  Z <- matrix(1, nrow = 2, ncol = 2)
  B <- matrix(1:6, nrow = 2, ncol = 3)
  P <- 3; K <- 2
  beta <- matrix(1:(P*K), nrow = P*K, ncol = 1)
  alpha_i <- list(rep(1, K), rep(1, K))

  partials_il <- get_partials_il(i, l, Y0, Z, B, beta, alpha_i, mi_vec, i_index, P, K)
  expect_type(partials_il, "double")
  expect_equal(dim(partials_il), c(P*K, K * mi_vec[i]))
})

test_that("get_partials_l_list returns list of correct matrices", {
  l <- 0
  mi_vec <- c(2, 3)
  i_index <- c(0, 2)
  Y0 <- 1:5
  Z <- matrix(1, nrow = 5, ncol = 2)
  B <- matrix(1:15, nrow = 5, ncol = 3)
  P <- 3; K <- 2
  beta <- matrix(1:(P*K), nrow = P*K, ncol = 1)
  alpha <- list(
    list(rep(1, K), rep(1, K)),
    list(rep(1, K), rep(1, K), rep(1, K))
  )

  partials_l <- get_partials_l_list(Y0, l, mi_vec, i_index, beta, alpha, Z, B, P, K)

  expect_length(partials_l, length(mi_vec))
  expect_equal(dim(partials_l[[1]]), c(P*K, K * mi_vec[1]))
  expect_equal(dim(partials_l[[2]]), c(P*K, K * mi_vec[2]))
})
