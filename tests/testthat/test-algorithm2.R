test_that("Algorithm 2", {
  sim <- base_sim()

  lp <- 1
  D <- diag(10)
  result2 <- algorithm2(Y = sim$Y,
                        Z = sim$Z,
                        is = sim$is,
                        mis = sim$mis,
                        lp = 1,
                        B = sim$B,
                        beta = sim$beta,
                        D = D,
                        gammas = rep(1,3),
                        phi = .5)
})
