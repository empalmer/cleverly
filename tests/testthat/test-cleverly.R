test_that("GEEwrapper works", {

  Y_mat <-  readRDS(test_path("test_data", "Y.rds"))
  sid <- rep(1:5, each = 3)
  tid <- rep(1:3, 5)
  Z_mat <- readRDS(test_path("test_data", "Z.rds"))

  # Test the wrapper function

  result <- GEEclevR(Y = Y_mat,
                     Z = Z_mat,
                     subject_ids = sid,
                     time = tid,
                     lp = 1,
                     gamma = rep(1, 3),
                     d = 3,
                     nknots = 3,
                     order = 3,
                     tol = 1e6,
                     smax = 2)

  expect_type(result, "list")
})
