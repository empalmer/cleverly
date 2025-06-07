test_that("Simulation Z", {
  skip("Skip - used as test file for cleverly")
  # Generate simulation data
  set.seed(127)
  sim <- simulation_data(n = 20,
                          range_start = 5000,
                          range_end = 20000,
                          nknots = 3,
                          K = 12,
                          order = 3,
                          user_var = 1000000,
                          cor_str = "CON-d",
                          rho = 0.95,
                          prob1 = .5,
                         baseline_fxns = list(
                           function(t) cos(2 * pi * t),
                           function(t) cos(2 * pi * t),
                           function(t) cos(2 * pi * t),
                           function(t) cos(2 * pi * t),
                           function(t) 2 * sin(pi * t) - 1,
                           function(t) 2 * sin(pi * t) - 1,
                           function(t) 2 * sin(pi * t) - 1,
                           function(t) 2 * sin(pi * t) - 1,
                           function(t) 2 - 2 * t,
                           function(t) 2 - 2 * t,
                           function(t) 2 - 2 * t,
                           function(t) 2 - 2 * t
                         ),
                          # Slope functions
                          slope_fxns = list(
                            function(t) 2 - t,
                            function(t) 2 * sin(pi * t),
                            function(t) .5,
                            function(t) t^2,
                            function(t)  -.5 ,
                            function(t) 2 * t - 2,
                            function(t) t,
                            function(t) 2 - t,
                            function(t) 1,
                            function(t) -.75,
                            function(t) .75,
                            function(t) 2 * t))

  # This is the one where it did not work in simulation
  # Setting 1: binary baseline
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/cleverly_simulation/simulation_data_results/sim_data/Z_binary_baseline/sim_data_47.rds")
  # Setting 2: binary slope
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/cleverly_simulation/simulation_data_results/sim_data/Z_binary_slope/sim_data_27.rds")
  # Setting 3: cont baseline
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/cleverly_simulation/simulation_data_results/sim_data/Z_cont_baseline/sim_data_1.rds")
  # Setting 4: cont slope
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/cleverly_simulation/simulation_data_results/sim_data/Z_cont_slope/sim_data_27.rds")

  # Visualize simulated data
  plot_sim_data(sim, Z_type = "binary")

  true_cluster <- rep(1:3, each = 4)

  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))
  Z <- sim$Z

  res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  time = time,
                  cluster_index = 1,
                  cor_str = "IND",
                  # Hyperparameters
                  gammas = c(1,1),
                  theta = 500,
                  parralel = F,
                  psi_min = 500,
                  psi_max = 2200,
                  npsi = 1,
                  # Iterations max
                  run_min = 3,
                  max_admm_iter = 100,
                  max_outer_iter = 3,
                  max_2_iter = 100,
  ) %>%
    get_cluster_diagnostics(true_cluster)


  res$y_hat_baseline


  plot_clusters(res,
                response_names = 1:12,
                curve_type = "slope",
                Y_counts = dplyr::select(Y, -c(time, individual)))

  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = 1:12,
                   curve_type = "slope",
                   Y_counts = dplyr::select(Y, -c(time, individual)))

  plot_clusters(res,
                response_names = 1:12,
                curve_type = "baseline",
                Y_counts = dplyr::select(Y, -c(time, individual)))

  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = 1:12,
                   curve_type = "baseline",
                   Y_counts = dplyr::select(Y, -c(time, individual)))



  plot_BIC(res, BIC_type = "BIC", psis = seq(600, 1400, length.out = 2))
  plot_BIC(res, BIC_type = "BIC_ra_group", psis = seq(600, 1400, length.out = 2))

  res$possible_clusters
  res$psi
  res$phi
  res$rho

  res$BIC

  res$clusters



  plot_initial_fit(res, K = 12)


})



test_that("CLR", {
  skip("Skip - used as test file for cleverly")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_binary_baseline/sim_data_27.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_cont_baseline/sim_data_27.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_binary_slope/sim_data_27.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_binary_baseline/sim_data_27.rds")


  true_cluster <- rep(1:3, each = 4)
  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z",
    "time",
    "individual"))
  Z <- sim$Z

  # Calculate B-spline basis based on time for each subject/time
  B <- get_B(time = sim$time,
             order = 3,
             nknots = 3)

  clr_res <- CLR_cluster(Y = Y,
                         Z = cbind(1, Z),
                         time = sim$time,
                         B = B,
                         lp = 0,
                         K = 12,
                         P = 6,
                         M = nrow(Y),
                         cluster_method = "gmm") %>%
    get_cluster_diagnostics(true_cluster)


  clr_res$cluster_diagnostics
  clr_res$clusters



  clr_res$y_hat_baseline

  plot_clusters(clr_res,
                response_names = 1:12,
                Z_type = "continuous",
                y_type = "y_hat_baseline",
                scales = "fixed")
  plot_sim_data(sim, Z_type = "continuous")


})

test_that("SMALL", {
  skip("Skip - used as test file for cleverly")
  # Generate simulation data
  set.seed(127)
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

  # Visualize simulated data
  plot_sim_data(sim)

  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))
  Z <- sim$Z

  res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  time = time,
                  lp = 0,
                  cor_str = "IND",
                  # Hyperparameters
                  gammas = c(1,1),
                  psi_min = 400,
                  npsi = 1,
                  # Iterations max
                  max_admm_iter = 10,
                  max_outer_iter = 5,
                  max_2_iter = 10,
  )


  res %>%
    get_cluster_diagnostics(true_cluster = rep(1:3, each = 4))
  res$clusters
  # Diagnostic plots:
  plot_clusters(res = res,
                response_names = LETTERS[1:12])

  res$phi
  res$rho


  plot_initial_fit(res, K = 12)

  # Get cluster membership
  cluster_df <- data.frame(
    K = factor(1:12, levels = 1:12),
    cluster = factor(res$clusters$membership))
  knitr::kable(table(cluster_df$cluster, rep(1:3, each = 4)))


  get_corR("AR1-d", mi = 3, K = 2, rho = .4)


})


