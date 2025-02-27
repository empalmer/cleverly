test_that("Update v", {
  sim <- base_sim()

  K <- 4
  Kappa <- t(combn(K,2))
  P <- 6
  Kappa_size <- nrow(Kappa)

  v <- update_v(beta = sim$beta,
                lp = 2,
                lambda = numeric(Kappa_size*P),
                Kappa = Kappa,
                P = P,
                gammas = rep(1,3),
                tau = 8/1000,
                theta = 300,
                psi = 1)

  expect_length(v, Kappa_size * P)
})

test_that("Update lambda", {
  sim <- base_sim()

  K <- 4
  Kappa <- t(combn(K,2))
  P <- 6
  lp <- 2
  Kappa_size <- nrow(Kappa)

  lambda <- numeric(Kappa_size*P)
  A <- get_A(Kappa = Kappa, K = K, P = P)

  Kappa_size <- nrow(Kappa)
  v <- update_v(beta = sim$beta,
                lp = lp,
                lambda = lambda,
                Kappa = Kappa,
                P = P,
                gammas = rep(1,3),
                tau = 8/1000,
                theta = 300,
                psi = 1)
  beta_lp <- sim$beta[, lp + 1]
  lambda <- update_lambda(beta_lp = beta_lp,
                          v = v,
                          lambda = lambda,
                          A = A,
                          theta = 1)

  expect_length(lambda, Kappa_size * P)
})


test_that("Update beta_lp", {
  sim <- base_sim()

  K <- 4
  Kappa <- t(combn(K,2))
  P <- 6
  lp <- 2
  Kappa_size <- nrow(Kappa)

  lambda <- numeric(Kappa_size*P)
  A <- get_A(Kappa = Kappa, K = K, P = P)
  B <- get_B(time = sim$time,
             order = 3,
             nknots = 3)

  Kappa_size <- nrow(Kappa)
  v <- update_v(beta = sim$beta,
                lp = lp,
                lambda = lambda,
                Kappa = Kappa,
                P = P,
                gammas = rep(1,3),
                tau = 8/1000,
                theta = 300,
                psi = 1)
  beta <- update_beta_admm(Y = sim$Y,
                           beta = sim$beta,
                           lp = lp,
                           v = v,
                           lambda = lambda,
                           theta = 1,
                           gammas = rep(1, 3),
                           D = get_D(K = K,
                                     d = 2,
                                     order = 3,
                                     nknots = 3),
                           A = A,
                           Z = sim$Z,
                           B = B,
                           phi = 1,
                           C = 10,
                           mi_vec = sim$mi_vec)

  expect_equal(dim(beta), dim(sim$beta))
  # Only the beta for lp should be different
  expect_equal(beta[,-(lp + 1)], sim$beta[,-(lp + 1)])
})


test_that("ADMM algorithm 3", {
  sim <- base_sim()

  K <- 4
  Kappa <- t(combn(K,2))
  P <- 6
  lp <- 2
  Kappa_size <- nrow(Kappa)

  A <- get_A(Kappa = Kappa, K = K, P = P)
  B <- get_B(time = sim$time,
             order = 3,
             nknots = 3)
  C <- 10
  D <- get_D(K = K,
             d = 2,
             order = 3,
             nknots = 3)

  beta_lp <- algorithm3(Y = sim$Y,
                        Z = sim$Z,
                        lp = lp,
                        B = B,
                        beta = sim$beta,
                        A = A,
                        P = P,
                        C = C,
                        D = D,
                        Kappa = Kappa,
                        mi_vec = sim$mi_vec,
                        gammas = rep(.001,3),
                        tau = 8/1000,
                        theta = 300,
                        psi = .001,
                        phi = 1,
                        max_admm_iter = 5)


  beta_lp$beta_admm_track[[1]][,3]
  beta_lp$beta_admm_track[[2]][,3]
  beta_lp$beta_admm_track[[4]][,3]



})

