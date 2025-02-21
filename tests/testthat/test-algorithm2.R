test_that("Algorithm 2", {
  sim <- base_sim()

  lp <- 2
  D <- get_D(K = 4,
             d = 2,
             order = 3,
             nknots = 3)
  beta_step <- algorithm2(Y = sim$Y,
                        Z = sim$Z,
                        is = sim$is,
                        mi_vec = sim$mi_vec,
                        lp = lp,
                        B = sim$B,
                        beta = sim$beta,
                        D = D,
                        gammas = rep(.01,3),
                        phi = .5,
                        C = 1)


})
