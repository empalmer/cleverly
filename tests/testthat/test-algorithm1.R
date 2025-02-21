test_that("Algorithm 1", {
  sim <- base_sim()

  result1 <- algorithm1(Y = sim$Y,
                        Z = sim$Z,
                        time = sim$time,
                        mi_vec = sim$mi_vec,
                        lp = 0,
                        gammas = rep(.01, 3),
                        psi = 1,
                        phi = 1,
                        tau = 8/1000,
                        theta = 300,
                        d = 2,
                        nknots = 3,
                        order = 3,
                        tol = 1e6,
                        max_outer_iter = 2,
                        max_admm_iter = 2)

  # Test beta is the expected dimension
  expect_equal(dim(result1$beta), c(6*4, 3))

})
