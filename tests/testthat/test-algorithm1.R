test_that("Algorithm 1", {
  sim <- simulation_data(n = 2,
                         range_start = 5000,
                         range_end = 6000,
                         nknots = 3,
                         K = 2,
                         order = 3,
                         maxt = .15,
                         user_var = 100,
                         cor_str = "CON-d",
                         rho = 0.5,
                         prob1 = .5,
                         baseline_fxns = list(
                           function(t) t,
                           function(t) -3 * t + 3
                         ),
                         slope_fxns = list(
                           function(t) 2,
                           function(t) t
                         ))

  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))
  Z <- sim$Z
  mi_vec <- get_mi_vec(Y_user, subject_ids, time)$mi

  result1 <- algorithm1(Y,
                        Z,
                        time,
                        mi_vec = mi_vec,
                        lp = 0,
                        cor_str = "IND",
                        gammas = c(1,1),
                        psi = 600,
                        tau = 0.005,
                        theta = 300,
                        C = 100,
                        d = 2,
                        nknots = 3,
                        order = 3,
                        epsilon_b = 1,
                        epsilon_r = 1,
                        epsilon_d = 1,
                        epsilon_2 = 1,
                        run_min = 1,
                        max_outer_iter = 1,
                        max_admm_iter = 1,
                        max_2_iter = 1)

  # Test beta is the expected dimension
  # expect_equal(dim(result1$beta), c(6*4, 3))

})
