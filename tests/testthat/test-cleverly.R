test_that("cleverly wrapper works", {

  skip("TODO")

  Y_mat <-  readRDS(test_path("test_data", "Y.rds"))
  sid <- rep(1:5, each = 3)
  tid <- rep(1:3, 5)
  Z_mat <- readRDS(test_path("test_data", "Z.rds"))
  Z <- Z_mat[,-1]

  # Test the wrapper function

  result <- cleverly(Y = Y_mat,
                     Z = Z,
                     subject_ids = sid,
                     time = tid,
                     lp = 1,
                     gammas = rep(1, 3),
                     psi = 1,
                     tau = 8/100,
                     theta = 1,
                     d = 2,
                     nknots = 3,
                     order = 3,
                     epsilon_b = 1e-6,
                     epsilon_r = 1e-6,
                     epsilon_d = 1e-6,
                     max_outer_iter = 2,
                     max_admm_iter = 2)

  expect_type(result, "list")
})
