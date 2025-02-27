test_that("Algorithm 1", {
  sim <- base_sim()

  result1 <- algorithm1(Y = sim$Y,
                        Z = sim$Z,
                        time = sim$time,
                        mi_vec = sim$mi_vec,
                        lp = 0,
                        gammas = rep(.001, 3),
                        psi = .01,
                        tau = 8/1000,
                        theta = 300,
                        C = 10,
                        d = 2,
                        nknots = 3,
                        order = 3,
                        epsilon_b = 1e-6,
                        epsilon_r = 1e-6,
                        epsilon_d = 1e-6,
                        max_outer_iter = 2,
                        max_admm_iter = 2)

  # Test beta is the expected dimension
  expect_equal(dim(result1$beta), c(6*4, 3))

})
