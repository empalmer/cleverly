

# Variance ----------------------------------------------------------------

test_that("get_V_ijj returns correct dimensions and values", {
  K <- 3
  alpha_ij <- rep(0.2, K)
  phi <- 2
  Y_ij0 <- 10

  result <- get_V_ijj(Y_ij0, phi, alpha_ij)

  expect_equal(dim(result), c(K, K))
  expect_true(is.matrix(result))
})


# Inverse variance --------------------------------------------------------

test_that("get_Vi_inv with IND", {
  i <- 1
  K <- 2
  mi_vec <- c(2)
  Y <- matrix(1:8, ncol = K)
  Y0 <- rowSums(Y)
  phi <- 1
  beta <- matrix(0.1, nrow = 2*K, ncol = 1)
  Z <- matrix(1, nrow = 2, ncol = 1)
  B <- matrix(1, nrow = 2, ncol = 2)
  i_index <- c(0, 2)
  alpha <- list(list(rep(0.2, K), rep(0.2, K)))

  result <- get_Vi_inv(i, Y, Y0, alpha, mi_vec, i_index, phi, beta, Z, B, K, "IND", rho_cor = NULL)

  expect_equal(dim(result), c(K * mi_vec[i], K * mi_vec[i]))
  expect_true(is.matrix(result))
})

test_that("get_Vi_inv with CON", {
  i <- 1
  K <- 2
  mi_vec <- c(2)
  Y <- matrix(1:8, ncol = K)
  Y0 <- rowSums(Y)
  phi <- 1
  beta <- matrix(0.1, nrow = 2*K, ncol = 1)
  Z <- matrix(1, nrow = 2, ncol = 1)
  B <- matrix(1, nrow = 2, ncol = 2)
  i_index <- c(0, 2)
  alpha <- list(list(rep(0.2, K), rep(0.2, K)))

  result <- get_Vi_inv(i,
                       Y,
                       Y0,
                       alpha,
                       mi_vec,
                       i_index,
                       phi, beta, Z, B, K, "CON",
                       rho_cor = .5)

  expect_equal(dim(result), c(K * mi_vec[i], K * mi_vec[i]))
  expect_true(is.matrix(result))
})

test_that("get_Vi_inv with CON-d", {
  i <- 1
  K <- 2
  mi_vec <- c(2)
  Y <- matrix(1:8, ncol = K)
  Y0 <- rowSums(Y)
  phi <- 1
  beta <- matrix(0.1, nrow = 2*K, ncol = 1)
  Z <- matrix(1, nrow = 2, ncol = 1)
  B <- matrix(1, nrow = 2, ncol = 2)
  i_index <- c(0, 2)
  alpha <- list(list(rep(0.2, K), rep(0.2, K)))

  result <- get_Vi_inv(i,
                       Y,
                       Y0,
                       alpha,
                       mi_vec,
                       i_index,
                       phi, beta, Z, B, K, "CON-d",
                       rho_cor = .5)

  expect_equal(dim(result), c(K * mi_vec[i], K * mi_vec[i]))
  expect_true(is.matrix(result))
})
test_that("get_Vi_inv with AR1", {
  i <- 1
  K <- 2
  mi_vec <- c(2)
  Y <- matrix(1:8, ncol = K)
  Y0 <- rowSums(Y)
  phi <- 1
  beta <- matrix(0.1, nrow = 2*K, ncol = 1)
  Z <- matrix(1, nrow = 2, ncol = 1)
  B <- matrix(1, nrow = 2, ncol = 2)
  i_index <- c(0, 2)
  alpha <- list(list(rep(0.2, K), rep(0.2, K)))

  result <- get_Vi_inv(i,
                       Y,
                       Y0,
                       alpha,
                       mi_vec,
                       i_index,
                       phi, beta, Z, B, K, "AR1",
                       rho_cor = .5)

  expect_equal(dim(result), c(K * mi_vec[i], K * mi_vec[i]))
  expect_true(is.matrix(result))
})

test_that("get_Vi_inv with AR1-d", {
  i <- 1
  K <- 2
  mi_vec <- c(2)
  Y <- matrix(1:8, ncol = K)
  Y0 <- rowSums(Y)
  phi <- 1
  beta <- matrix(0.1, nrow = 2*K, ncol = 1)
  Z <- matrix(1, nrow = 2, ncol = 1)
  B <- matrix(1, nrow = 2, ncol = 2)
  i_index <- c(0, 2)
  alpha <- list(list(rep(0.2, K), rep(0.2, K)))

  result <- get_Vi_inv(i,
                       Y,
                       Y0,
                       alpha,
                       mi_vec,
                       i_index,
                       phi, beta, Z, B, K, "AR1-d",
                       rho_cor = .5)

  expect_equal(dim(result), c(K * mi_vec[i], K * mi_vec[i]))
  expect_true(is.matrix(result))
})

test_that("get_V_inv returns list of correct size and matrices", {
  K <- 2
  mi_vec <- c(2, 1)
  n <- length(mi_vec)
  M <- sum(mi_vec)
  Y <- matrix(1:(M * K), ncol = K)
  Y0 <- rowSums(Y)
  phi <- 1
  beta <- matrix(0.1, nrow = 2*K, ncol = 1)
  Z <- matrix(1, nrow = M, ncol = 1)
  B <- matrix(1, nrow = M, ncol = 2)
  i_index <- c(0, 2)
  alpha <- list(
    list(rep(0.2, K), rep(0.2, K)),
    list(rep(0.2, K))
  )

  result <- get_V_inv(Y, Y0, alpha, mi_vec, i_index, phi, beta, Z, B, K, "IND", rho_cor = NULL)

  expect_equal(length(result), n)
  expect_true(all(sapply(result, is.matrix)))
})



# Correlation structures --------------------------------------------------


test_that("get_corR works for IND structure", {
  result <- get_corR("IND", mi = 2, K = 2, rho = 0.5)
  expect_equal(result, list())
})

test_that("get_corR works for CON structure", {
  mi <- 3; K <- 2; rho <- 0.5
  result <- get_corR("CON", mi, K, rho)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(mi * K, mi * K))
  expect_true(all(diag(result) == 0))
})

test_that("get_corR works for AR1 structure", {
  mi <- 4; K <- 1; rho <- 0.3
  result <- get_corR("AR1", mi, K, rho)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(mi * K, mi * K))
  expect_equal(result[1, 2], rho)
  expect_equal(result[1, 3], rho^2)
  expect_equal(result[1, 1], 0)  # diagonal should be 0
})

test_that("get_corR works for CON-d structure", {
  mi <- 3; K <- 2; rho <- 0.5
  result <- get_corR("CON-d", mi, K, rho)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(mi * K, mi * K))
  expect_equal(unique(diag(result)), 0)
  expect_true(all(result[result != 0] == rho))
})

test_that("get_corR works for AR1-d structure", {
  mi <- 3; K <- 2; rho <- 0.2
  result <- get_corR("AR1-d", mi, K, rho)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(mi * K, mi * K))
  expect_equal(result[1, K + 1], rho)
  expect_equal(result[1, 1 + 2 * K], rho^2)
})

test_that("get_rho returns 0 for IND", {
  pearson_residuals <- list(matrix(rnorm(4), ncol = 1))
  phi <- 1
  K <- 1
  mi_vec <- c(2)
  M <- 2
  cor_str <- "IND"
  cor_blocks <- list()  # unused
  j1_j2_list <- list()

  result <- get_rho(pearson_residuals, phi, K, mi_vec, M, cor_str, cor_blocks, j1_j2_list)
  expect_equal(result, 0)
})

test_that("get_rho returns numeric within bounds for CON", {
  K <- 1
  mi_vec <- c(3)
  M <- sum(mi_vec)
  phi <- 1

  pearson_residuals <- list(rnorm(K * mi_vec[1]))
  cor_blocks <- list(get_corR("CON", mi_vec[1], K, 1))
  j1_j2_list <- list()

  result <- get_rho(pearson_residuals, phi, K, mi_vec, M, "CON", cor_blocks, j1_j2_list)
  expect_true(is.numeric(result))
  expect_true(result >= -1 && result <= 1)
})

test_that("get_rho returns numeric within bounds for AR1", {
  mi_vec <- c(4)
  K <- 1
  M <- sum(mi_vec)
  phi <- 1

  residuals_i <- rnorm(mi_vec[1] * K)
  pearson_residuals <- list(residuals_i)

  cor_blocks <- list(get_corR("AR1", mi_vec[1], K, 1))

  j1s <- matrix(rep(1:mi_vec[1], each = mi_vec[1]), nrow = mi_vec[1])
  j2s <- matrix(rep(1:mi_vec[1], each = mi_vec[1]), ncol = mi_vec[1], byrow = TRUE)
  j1_j2_list <- list(j1s - j2s)

  result <- get_rho(pearson_residuals, phi, K, mi_vec, M, "AR1", cor_blocks, j1_j2_list)
  expect_true(is.numeric(result))
  expect_true(result >= -1 && result <= 1)
})

test_that("get_rho handles CON-d structure", {
  mi_vec <- c(3)
  K <- 2
  M <- sum(mi_vec)
  phi <- 1

  residuals_i <- rnorm(mi_vec[1] * K)
  pearson_residuals <- list(residuals_i)
  cor_blocks <- list(get_corR("CON-d", mi_vec[1], K, 1))
  j1_j2_list <- list()

  result <- get_rho(pearson_residuals, phi, K, mi_vec, M, "CON-d", cor_blocks, j1_j2_list)
  expect_true(is.numeric(result))
  expect_true(result >= -1 && result <= 1)
})
