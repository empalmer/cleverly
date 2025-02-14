

# Various -----------------------------------------------------------------

test_that("extract yij", {
  Y <- readRDS(test_path("test_data", "Y.rds"))
  expect_equal(get_Y_ij0(i = 1, j = 1, Y = Y, mi = 3), 100)
})


test_that("B", {
  B <- readRDS(test_path("test_data", "B.rds"))
  expect_length(get_B_ij(i = 1, j = 1, B = B, mi = 3), 6)

})

test_that("Z", {
  Z <- readRDS(test_path("test_data", "Z.rds"))
  expect_equal(get_Z_ijl(i = 1, j = 2, l = 2, Z = Z), 1)
})


test_that("Z0 all 1", {
  Z <- readRDS(test_path("test_data", "Z.rds"))
  # tests for the first time point of each sample
  for (i in 1:nrow(Z)) {
    expect_equal(get_Z_ijl(i = i, j = 1, l = 0, Z = Z), 1)
  }
})




# Alpha -------------------------------------------------------------------

test_that("Check alpha_ijk", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4

  alpha_ijk <- get_alpha_ijk(i = 1, j = 1, k = 1, beta = beta, Z = Z, B = B, mi = 3)
  # Check dimensions
  # Should be one value
  expect_length(alpha_ijk, 1)
})


test_that("Check alpha_ij", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data","Z.rds"))
  B <- readRDS(test_path("test_data","B.rds"))
  K <- 4

  alpha_ij <- get_alpha_ij(i = 1, j = 1, beta = beta, Z = Z, B = B, K = K, mi = 3)
  # Should be of dimension K
  expect_length(alpha_ij, K)

  expect_equal(alpha_ij[1], 1.204884947)
})



# Mus ---------------------------------------------------------------------

test_that("Check mu_ij", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4

  # Test if alpha is supplied
  mu_ij_alpha <- get_mu_ij(Y_ij0 = 100,
                     alpha_ij = get_alpha_ij(i = 1, j = 1, beta = beta, Z = Z, B = B, K = K, mi = 3))
  # Should be of dimension K
  expect_length(mu_ij_alpha, K)
  expect_equal(round(mu_ij_alpha[1]), 30.0)

  # Test if alpha is not supplied
  mu_ij <- get_mu_ij(Y_ij0 = 100, i = 1, j = 1, beta = beta, Z = Z, B = B, K = K, mi = 3)
  expect_length(mu_ij, K)
  expect_equal(round(mu_ij[1]), 30.0)
})


test_that("Check mu_i", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3

  # Test if alpha is supplied
  mu_i <- get_mu_i(i = 1, mi = 3, Y = Y, beta = beta, Z = Z, B = B, K = K)

  # Should be of dimension K
  expect_length(mu_i, K*mi)
})



# Variance -----------------------------------------------------

test_that("Check Uij", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4

  U_ij <- get_U_ij(Y_ij0 = 100,
                   alpha_ij = get_alpha_ij(i = 1,
                                           j = 2,
                                           beta = beta,
                                           Z = Z,
                                           B = B,
                                           K = K,
                                           mi = 3))
  # Check is square
  expect_equal(dim(U_ij)[1], dim(U_ij)[2])

  # Check right dimension:
  expect_equal(dim(U_ij)[1], K)
})





# Check case where L = 0 --------------------------------------------------


