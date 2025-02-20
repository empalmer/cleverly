test_that("Algorithm 1", {
  set.seed(124)
  time_list <- sim_timepoints(n = 5)
  mis <- time_list$mis
  time <- time_list$X$time
  subject_ids <- rep(1:5, mis)

  K <- 4
  Z <- sim_Z(mis)

  B <- get_B(time, order = 3, nknots = 3)

  # 2 clusters
  betaC1 <- matrix(c(rep(c(1, 1, 1, 1, 1, 1), 2),     #l = 0
                     rep(c(-2, -2, -2, -2, -2, -2),2),  #l = 1
                     rep(c(1, 2, 3, 4, 5, 6), 2)), ncol = 3)#l = 2
  betaC2 <- .5*betaC1
  beta <- rbind(betaC1, betaC2)


  Y <- sim_Y(beta = beta,
             Z = Z,
             B = B,
             K = K,
             mis = mis)

  result1 <- algorithm1(Y = Y,
                        Z = Z,
                        is = is,
                        time = time,
                        mis = mis,
                        lp = 1,
                        gamma = rep(1, 3),
                        phi = .5,
                        d = 3,
                        nknots = 3,
                        order = 3,
                        tol = 1e6,
                        smax = 2)



})
