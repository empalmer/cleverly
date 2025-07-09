test_that("Case 1: Binary baseline", {
  # Generate simulation data
  set.seed(127)
  sim <- simulation_data(n = 20,
                         range_start = 5000,
                         range_end = 20000,
                         nknots = 3,
                         K = 12,
                         order = 3,
                         user_var = 1000,
                         cor_str = "CON-d",
                         rho = 0.5,
                         prob1 = .5,
                         Z_type = "binary",
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

  # Visualize simulated data
  plot_sim_data(sim)  %>%
    expect_s3_class("ggplot")

  true_cluster <- rep(1:3, each = 4)

  Y <- sim$Y
  Z <- sim$Z

  # Testing with only 1 psi so we should get non-unique clusters warning
  expect_warning(
    res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  time = time,
                  cluster_index = 0,
                  cor_str = "IND",
                  theta = 500,
                  parralel = F,
                  psi_min = 500,
                  psi_max = 501,
                  npsi = 1,
                  # Iterations max
                  max_admm_iter = 2,
                  max_2_iter = 2) %>%
    get_cluster_diagnostics(true_cluster))


  res$clusters

  expect_null(res$error)

  plot_initial_fit(res,
                   response_names = LETTERS[1:12]) %>%
    expect_s3_class("ggplot")

  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "baseline") %>%
    expect_s3_class("ggplot")
  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "refit_baseline") %>%
    expect_s3_class("ggplot")
  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "baseline") %>%
    expect_s3_class("ggplot")

})


test_that("Case 2: Binary slope", {
  # Generate simulation data
  set.seed(127)
  sim <- simulation_data(n = 20,
                         range_start = 5000,
                         range_end = 20000,
                         nknots = 3,
                         K = 12,
                         order = 3,
                         user_var = 10000,
                         cor_str = "CON-d",
                         rho = 0.5,
                         prob1 = .5,
                         Z_type = "binary",
                         slope_fxns = list(
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
                         baseline_fxns = list(
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


  # Visualize simulated data
  plot_sim_data(sim)  %>%
    expect_s3_class("ggplot")

  true_cluster <- rep(1:3, each = 4)

  Y <- sim$Y
  Z <- sim$Z

  # Testing with only 1 psi so we should get non-unique clusters warning
  expect_warning(res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  time = time,
                  cluster_index = 1,
                  cor_str = "IND",
                  theta = 500,
                  parralel = F,
                  psi_min = 500,
                  psi_max = 501,
                  npsi = 1,
                  # Iterations max
                  max_admm_iter = 2,
                  max_2_iter = 2) %>%
    get_cluster_diagnostics(true_cluster))

  res$clusters

  plot_initial_fit(res,
                   response_names = LETTERS[1:12])  %>%
    expect_s3_class("ggplot")

  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "slope")  %>%
    expect_s3_class("ggplot")

  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "slope")  %>%
    expect_s3_class("ggplot")

})


test_that("Case 3: Continuous baseline", {
  # Generate simulation data
  set.seed(127)
  sim <- simulation_data(n = 20,
                         range_start = 5000,
                         range_end = 20000,
                         nknots = 3,
                         K = 12,
                         order = 3,
                         user_var = 1000,
                         cor_str = "CON-d",
                         rho = 0.5,
                         prob1 = .5,
                         Z_type = "continuous",
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


  # Visualize simulated data
  plot_sim_data(sim)

  true_cluster <- rep(1:3, each = 4)

  Y <- sim$Y
  Z <- sim$Z

  # Testing with only 1 psi so we should get non-unique clusters warning
  expect_warning(
    res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  time = time,
                  cluster_index = 0,
                  cor_str = "IND",
                  theta = 500,
                  parralel = F,
                  psi_min = 500,
                  psi_max = 501,
                  npsi = 1,
                  # Iterations max
                  max_admm_iter = 2,
                  max_2_iter = 2) %>%
    get_cluster_diagnostics(true_cluster))

  res$clusters

  plot_initial_fit(res,
                   response_names = LETTERS[1:12],
                   Z_name = "Z")  %>%
    expect_s3_class("ggplot")

  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "baseline")  %>%
    expect_s3_class("ggplot")
  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "baseline")  %>%
    expect_s3_class("ggplot")


})


test_that("Case 4: Continuous slope", {
  # Generate simulation data
  set.seed(127)
  sim <- simulation_data(n = 20,
                         range_start = 5000,
                         range_end = 20000,
                         nknots = 3,
                         K = 12,
                         order = 3,
                         user_var = 10000,
                         cor_str = "CON-d",
                         rho = 0.5,
                         prob1 = .5,
                         Z_type = "continuous",
                         slope_fxns = list(
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
                         baseline_fxns = list(
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


  # Visualize simulated data
  plot_sim_data(sim)

  true_cluster <- rep(1:3, each = 4)

  Y <- sim$Y
  Z <- sim$Z

  # Testing with only 1 psi so we should get non-unique clusters warning
  expect_warning(
    res <- cleverly(Y = Y,
                    Z = Z,
                    subject_ids = individual,
                    time = time,
                    cluster_index = 1,
                    cor_str = "IND",
                    theta = 500,
                    parralel = F,
                    psi_min = 500,
                    psi_max = 501,
                    npsi = 1,
                    # Iterations max
                    max_admm_iter = 2,
                    max_2_iter = 2) %>%
      get_cluster_diagnostics(true_cluster))

  res$clusters

  plot_initial_fit(res,
                   response_names = LETTERS[1:12]) %>%
    expect_s3_class("ggplot")

  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "slope",
                Y_counts = dplyr::select(Y, -c(time, individual))) %>%
    expect_s3_class("ggplot")

  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "slope",
                   Y_counts = dplyr::select(Y, -c(time, individual))) %>%
    expect_s3_class("ggplot")


})


test_that("Case 5: 2 Zs (BB)", {
  # Generate simulation data
  set.seed(127)
  sim <- simulation_data(n = 20,
                         range_start = 5000,
                         range_end = 20000,
                         nknots = 3,
                         K = 12,
                         order = 3,
                         user_var = 10000,
                         cor_str = "CON-d",
                         rho = 0.5,
                         prob1 = .5,
                         Z_type = "binary",
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


  # Visualize simulated data
  plot_sim_data(sim)

  true_cluster <- rep(1:3, each = 4)

  Y <- sim$Y
  Z <- sim$Z

  # Test 2 Zs
  Z <- cbind(Z, rnorm(nrow(Y), mean = 0, sd = 1)) # Add a random covariate
  colnames(Z) <- c("Z1", "Z2")

  # Testing with only 1 psi so we should get non-unique clusters warning
  expect_warning(
  res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  time = time,
                  cluster_index = 0,
                  cor_str = "IND",
                  theta = 500,
                  parralel = F,
                  psi_min = 500,
                  psi_max = 501,
                  npsi = 1,
                  # Iterations max
                  max_admm_iter = 2,
                  max_2_iter = 2) %>%
    get_cluster_diagnostics(true_cluster))

  res$clusters

  plot_initial_fit(res,
                   response_names = LETTERS[1:12],
                   Z_name = "Z1") %>%
    expect_s3_class("ggplot")

  # Doesnt work well when Z2 is continuous, but still outputs plot with right shading.
  plot_initial_fit(res,
                   response_names = LETTERS[1:12],
                   Z_name = "Z2") %>%
    expect_s3_class("ggplot")


  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "baseline",
                Z_name = "Z1") %>%
    expect_s3_class("ggplot")
  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "baseline",
                Z_name = "Z2") %>%
    expect_s3_class("ggplot")

  # Binary Z:
  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "baseline",
                   Z_name = "Z1") %>%
    expect_s3_class("ggplot")
  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "slope",
                   Z_name = "Z1") %>%
    expect_s3_class("ggplot")

  # Continuous Z
  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "baseline",
                   Z_name = "Z2") %>%
    expect_s3_class("ggplot")
  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "slope",
                   Y_counts = dplyr::select(Y, -c(time, individual)),
                   Z_name = "Z2") %>%
    expect_s3_class("ggplot")
})



test_that("Case 5: Categorical Zs (BB)", {
  # Generate simulation data
  set.seed(127)
  sim <- simulation_data(n = 20,
                         range_start = 5000,
                         range_end = 20000,
                         nknots = 3,
                         K = 12,
                         order = 3,
                         user_var = 10000,
                         cor_str = "CON-d",
                         rho = 0.5,
                         prob1 = .5,
                         Z_type = "binary",
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


  # Visualize simulated data
  plot_sim_data(sim)


  true_cluster <- rep(1:3, each = 4)

  Y <- sim$Y
  Z <- sim$Z
  Z <- ifelse(Z == 0, "bar", "foo")
  #Z <- factor(Z)

  expect_warning(
  res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  time = time,
                  cluster_index = 0,
                  cor_str = "IND",
                  theta = 500,
                  parralel = F,
                  psi_min = 500,
                  psi_max = 501,
                  npsi = 1,
                  # Iterations max
                  max_admm_iter = 2,
                  max_2_iter = 2) %>%
    get_cluster_diagnostics(true_cluster))

  res$clusters
  res$Z

  plot_initial_fit(res,
                   response_names = LETTERS[1:12],
                   Z_name = "Zfoo") %>%
    expect_s3_class("ggplot")

  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "baseline",
                Z_name = "Zfoo") %>%
    expect_s3_class("ggplot")

  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
                   curve_type = "slope",
                   Z_name = "Zfoo") %>%
    expect_s3_class("ggplot")
})


test_that("CLR", {
  skip("Skip - used as test file")
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

})

