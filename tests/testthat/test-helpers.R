

# Various -----------------------------------------------------------------

test_that("Y ij0", {
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))
  Y <- readRDS(test_path("test_data", "Y.rds"))
  Y_ij0 <- get_Y_ij0(i = 1, j = 1, Y = Y, i_index = i_index)
  expect_equal(Y_ij0, 100)
})


test_that("Y i", {
  Y <- readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)
  # Check dimensions
  Y_i_mat <- get_Y_i_mat(i = 1, mi_vec = mi_vec, Y = Y)
  expect_equal(dim(Y_i_mat)[1], mi)
  expect_equal(dim(Y_i_mat)[2], K)
})

test_that("Y i SIM ", {
  sim <- base_sim()

  i <- 3
  Y_i_mat <- get_Y_i_mat(i = i, mi_vec = sim$mi_vec, Y = sim$Y)
  Y_i_vec <- get_Y_i_vec(i = i, mi_vec = sim$mi_vec, Y = sim$Y)


  expect_equal(dim(Y_i_mat)[1], sim$mi_vec[i])
  expect_equal(dim(Y_i_mat)[2], sim$K)
})

test_that("Y vec and matrix", {
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)

  Y_i_mat <- get_Y_i_mat(i = 1, mi_vec = mi_vec, Y = Y)
  Y_i_vec <- get_Y_i_vec(i = 1, mi_vec = mi_vec, Y = Y)

  #j=1
  expect_equal(as.numeric(Y_i_mat[1, ]),
               as.numeric(Y_i_vec[1:K]))
  #j = 2
  expect_equal(as.numeric(Y_i_mat[2, ]),
               as.numeric(Y_i_vec[(K+1):(2*K)]))

})



test_that("Y and mu indexing match", {
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  Y_i <- get_Y_i_vec(i = 1,
                     mi_vec = mi_vec,
                     Y = Y)
  mu_i <- get_mu_i(i = 1,
                   mi_vec = mi_vec,
                   i_index = i_index,
                   Y = Y,
                   beta = beta,
                   Z = Z,
                   B = B,
                   K = K)

  expect_equal(length(Y_i), length(mu_i))
})




test_that("Test Y wrapper", {

  Y <-  readRDS(test_path("test_data", "Y.rds"))
  sid <- rep(1:5, each = 3)
  tid <- rep(1:3, 5)


  # Case where Y is a matrix and subject and time are external columns
  Y_mat <- get_Y_wrapper(
    Y = Y,
    subject_ids = sid,
    time_ids = tid)

  # Case where Y is a matrix and subject and time are external character columns.
  Y_mat_character <- get_Y_wrapper(
    Y = Y,
    subject_ids = rep(LETTERS[1:5], each = 3),
    time_ids = rep(letters[1:3], 5)
  )

  # Case where Y is a data frame and subject and time are exteral columns (character)
  Y_df_character <- get_Y_wrapper(
    Y = data.frame(Y),
    subject_ids = rep(LETTERS[1:5], each = 3),
    time_ids = rep(letters[1:3], 5)
  )


  # Case where Y is a data frame and subject and time are columns in Y
  Y_df_column <- get_Y_wrapper(
    Y = data.frame(Y, sid = sid, tid = tid),
    subject_ids = sid,
    time_ids = tid
  )


  expect_equal(Y_mat, Y_df_column)
  expect_equal(Y_mat_character, Y_df_character)


  # expect errors:
  expect_error(
    get_Y_wrapper(
      Y = data.frame(Y, sid = sid, tid = tid),
      subject_ids = sid,
      time_ids = 1:100
    ), "Invalid ID")

  expect_error(
    get_Y_wrapper(
      Y = 1:100,
      subject_ids = sid,
      time_ids = nope
    ), "Y must be either a data frame or a matrix")

})



test_that("Check mi_vec", {
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  sid <- rep(letters[1:5], each = 3)
  tid <- rep(1:3, 5)

  mi_vec <- get_mi_vec(Y, subject_ids = sid, time_ids = tid)
  expect_equal(mi_vec$mi, rep(3, 5))
  expect_equal(mi_vec$subject_id, letters[1:5])
})



test_that("B", {
  mi_vec <- rep(3, 5)
  B <- readRDS(test_path("test_data", "B.rds"))
  i_index <- c(0, cumsum(mi_vec))

  B_ij <- get_B_ij(i = 3, j = 2, B = B, i_index = i_index)

  expect_length(B_ij, 6)

})




test_that("Z", {
  Z <- readRDS(test_path("test_data", "Z.rds"))
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))
  Z_ijl <- get_Z_ijl(i = 5, j = 2,
                     l = 2, Z = Z,
                     i_index = i_index)
  expect_equal(Z_ijl, 1)
})



test_that("Z0 all 1", {
  Z <- readRDS(test_path("test_data", "Z.rds"))
  mi_vec <- rep(3, 5)
  n <- length(mi_vec)
  i_index <- c(0, cumsum(mi_vec))
  # tests for the first time point of each sample
  for (i in 1:n) {
    Z_ijl <- get_Z_ijl(i = i, j = 1,
                       l = 0, Z = Z,
                       i_index = i_index)
    expect_equal(Z_ijl, 1)
  }
})




# Alpha -------------------------------------------------------------------

test_that("Check alpha_ijk", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  alpha_ijk <- get_alpha_ijk(i = 1,
                             j = 1, k = 1,
                             beta = beta,
                             Z = Z, B = B,
                             i_index = i_index)

  # Check dimensions
  # Should be one value
  expect_length(alpha_ijk, 1)
})


test_that("Check alpha_ij", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data","Z.rds"))
  B <- readRDS(test_path("test_data","B.rds"))
  K <- 4
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))
  alpha_ij <- get_alpha_ij(i = 1, j = 1,
                           beta = beta,
                           Z = Z,
                           B = B,
                           K = K,
                           i_index = i_index)

  alpha_ij2 <- get_alpha_ij2(i = 1,
                           beta = beta,
                           Z = Z,
                           B = B,
                           K = K,
                           i_index = i_index,
                           mi = 3)
  # Should be of dimension K
  expect_length(alpha_ij, K)

  expect_equal(as.numeric(alpha_ij)[1], 1.204884947)
})

test_that("Check alpha list = alphaijk", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data","Z.rds"))
  B <- readRDS(test_path("test_data","B.rds"))
  K <- 4
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  i <- 2
  j <- 3
  k <- 2
  alpha <- get_alpha_list(beta, Z, B, K, i_index, mi_vec)
  alpha_ijk <- get_alpha_ijk(i = i,
                             j = j,
                             k, beta, Z, B, i_index)
  expect_equal(alpha_ijk, alpha[[i]][[j]][k])
})



# Mus ---------------------------------------------------------------------

test_that("Check mu_ij", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  # Test if alpha is supplied
  mu_ij_alpha <- get_mu_ij(Y_ij0 = 100,
                     alpha_ij = get_alpha_ij(i = 1,
                                             j = 1,
                                             beta = beta,
                                             Z = Z,
                                             B = B,
                                             K = K,
                                             i_index = i_index))
  # Should be of dimension K
  expect_length(mu_ij_alpha, K)
  expect_equal(round(as.numeric(mu_ij_alpha[1])), 30.0)

  # Test if alpha is not supplied
  mu_ij <- get_mu_ij(Y_ij0 = 100,
                     i = 1,
                     j = 1,
                     beta = beta,
                     Z = Z,
                     B = B,
                     K = K,
                     i_index = i_index)
  expect_length(mu_ij, K)
  expect_equal(round(as.numeric(mu_ij[1])), 30.0)
})


test_that("Check mu_i", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)

  # Test if alpha is supplied
  mu_i <- get_mu_i(i = 1,
                   mi_vec = mi_vec,
                   i_index = c(0, cumsum(mi_vec)),
                   Y = Y,
                   beta = beta,
                   Z = Z,
                   B = B,
                   K = K)

  # Should be of dimension K
  expect_length(mu_i, K*mi)
})


test_that("Check mu_i SIM", {
  sim <- base_sim()

  # Test if alpha is supplied
  i <- 3
  i_index <- c(0, cumsum(sim$mi_vec))
  mu_i <- get_mu_i(i = i,
                   mi_vec = sim$mi_vec,
                   i_index = i_index,
                   Y = sim$Y,
                   beta = sim$beta,
                   Z = sim$Z,
                   B = sim$B,
                   K = sim$K)

  # Should be of dimension K
  expect_length(mu_i, sim$K*sim$mi_vec[i])
})


# Variance -----------------------------------------------------

test_that("Check Uij", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  U_ij <- get_U_ij(alpha_ij = get_alpha_ij(i = 1,
                                           j = 2,
                                           beta = beta,
                                           Z = Z,
                                           B = B,
                                           K = K,
                                           i_index = i_index))
  # Check is square
  expect_equal(dim(U_ij)[1], dim(U_ij)[2])

  # Check right dimension:
  expect_equal(dim(U_ij)[1], K)
})


test_that("Check Vijj dimension", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  K <- 4
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  V_ijj <- get_V_ijj(Y_ij0 = 100,
                     phi = 1,
                     alpha_ij = get_alpha_ij(i = 1,
                                             j = 2,
                                             beta = beta,
                                             Z = Z,
                                             B = B,
                                             K = K,
                                             i_index = i_index))
  # Check is square
  expect_equal(dim( V_ijj)[1], dim( V_ijj)[2])

  # Check right dimension:
  expect_equal(dim(V_ijj)[1], K)
})

test_that("Check Vi dimension", {
  beta <- readRDS(test_path("test_data", "beta.rds"))
  Z <- readRDS(test_path("test_data", "Z.rds"))
  B <- readRDS(test_path("test_data", "B.rds"))
  Y <-  readRDS(test_path("test_data", "Y.rds"))
  K <- 4
  mi <- 3
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  V_i <- get_V_i(i = 1,
                 Y = Y,
                 phi = 1,
                 beta = beta,
                 Z = Z,
                 B = B,
                 K = K,
                 mi_vec = mi_vec,
                 i_index = i_index)
  # Check is square
  expect_equal(dim(V_i)[1], dim(V_i)[2])

  # Check right dimension:
  expect_equal(dim(V_i)[1], K * mi)

})


# Check case where L = 0 --------------------------------------------------


test_that("Check alpha_ijk l0", {
  beta <- readRDS(test_path("test_data", "betal0.rds"))
  Z <- readRDS(test_path("test_data", "Zl0.rds"))
  B <- readRDS(test_path("test_data", "Bl0.rds"))
  K <- 4
  mi_vec <- rep(3, 5)
  i_index <- c(0, cumsum(mi_vec))

  alpha_ijk <- get_alpha_ijk(i = 1,
                             j = 1,
                             k = 1,
                             beta = beta,
                             Z = Z,
                             B = B,

                             i_index = i_index)
  # Check dimensions
  # Should be one value
  expect_length(alpha_ijk, 1)
})





# Check cases when mi not all the same  -----------------------------------


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
