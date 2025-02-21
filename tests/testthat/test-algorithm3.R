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

  expect_length(beta, K * P)
})
