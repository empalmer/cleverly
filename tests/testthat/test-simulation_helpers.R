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
  plot_sim_data(sim)

  true_cluster <- rep(1:3, each = 4)

  Y <- sim$Y
  Z <- sim$Z

  # Test 2 Zs
  Z <- cbind(Z, rnorm(nrow(Y), mean = 0, sd = 1)) # Add a random covariate
  colnames(Z) <- c("Z1", "Z2")

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
                  max_admm_iter = 50,
                  max_2_iter = 50) %>%
    get_cluster_diagnostics(true_cluster)

  res$clusters


  res$y_hat_baseline

  plot_initial_fit(res,
                   response_names = LETTERS[1:12],
                   Z_col = "Z")


  plot_clusters(res,
                response_names = LETTERS[1:12],
                curve_type = "slope",
                Y_counts = dplyr::select(Y, -c(time, individual)))


  plot_cluster_differences(res, res)


  plot_one_cluster(res,
                   cluster_val = 1,
                   response_names = LETTERS[1:12],
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




  # Define slope and baseline functions
  slope_fxns <- list(
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
  )

  baseline_fxns <- list(
    function(t) 2 - t,
    function(t) 2 * sin(pi * t),
    function(t) 0.5,
    function(t) t^2,
    function(t) -0.5,
    function(t) 2 * t - 2,
    function(t) t,
    function(t) 2 - t,
    function(t) 1,
    function(t) -0.75,
    function(t) 0.75,
    function(t) 2 * t
  )

  # Create evaluation grid
  t_vals <- seq(0, 1, length.out = 100)

  # Build a long-format dataframe for ggplot
  plot_data <- map2_dfr(baseline_fxns, slope_fxns, ~{
    tibble(t = t_vals,
           baseline = .x(t_vals),
           slope = .y(t_vals),
           sum = .x(t_vals) + .y(t_vals))
  }, .id = "facet_id") %>%
    pivot_longer(cols = c("baseline", "sum"),
                 names_to = "curve_type",
                 values_to = "value") %>%
    mutate(facet_id = factor(paste("Function", facet_id),
                             levels = paste("Function", seq_along(baseline_fxns))),
           curve_type = recode(curve_type,
                               baseline = "Baseline",
                               sum = "Baseline + Slope"))

  # Plot
  ggplot(plot_data, aes(x = t, y = value, color = curve_type)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ facet_id, scales = "free_y") +
    labs(x = "t", y = "Value", color = "Curve Type",
         title = "Slope Curves f(t) and g(t) + f(t)") +
    theme(legend.position = "bottom")


  # Alpha plots:

  ggplot(plot_data, aes(x = t, y = exp(value), color = curve_type)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ facet_id, scales = "fixed") +
    labs(x = "t", y = "Value", color = "Curve Type",
         title = "Alpha values (exp(f(t) + Z g(t))) when Z = 0 and Z = 1") +
    theme(legend.position = "bottom")


  plot_data <- plot_data %>%
    mutate(alpha = exp(value)) %>%
    group_by(t, curve_type) %>%
    mutate(alpha0 = sum(alpha),
           alpha_ra = alpha / sum(alpha)) %>%
    ungroup()

  ggplot(plot_data, aes(x = t, y = alpha_ra, color = curve_type)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ facet_id, scales = "fixed") +
    labs(x = "t", y = "Value", color = "Curve Type",
         title = "Alpha/alpha0") +
    theme(legend.position = "bottom")

  ggplot(plot_data, aes(x = t, y = alpha/alpha0, color = curve_type)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ facet_id, scales = "fixed") +
    labs(x = "t", y = "Value", color = "Curve Type",
         title = "yhat ") +
    theme(legend.position = "bottom")


  plot_data <- plot_data %>%
    group_by(t, curve_type) %>%
    mutate(mu = alpha/alpha0,
           mu0 = mu/sum(mu)) %>%
    ungroup()


  ggplot(plot_data, aes(x = t, y = mu0, color = curve_type)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ facet_id, scales = "fixed") +
    labs(x = "t", y = "Value", color = "Curve Type",
         title = "yhat RA ") +
    theme(legend.position = "bottom")



  diff_data <- plot_data %>%
    dplyr::select(facet_id, t, mu0, curve_type) %>%
    pivot_wider(names_from = curve_type, values_from = mu0) %>%
    mutate(diff = `Baseline + Slope` - Baseline)

  ggplot(diff_data, aes(x = t, y = diff)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ facet_id, scales = "fixed") +
    labs(x = "t", y = "Value",
         title = "Diff  ")

  ggplot(diff_data, aes(x = t, y = exp(diff))) +
    geom_line(linewidth = 1) +
    facet_wrap(~ facet_id, scales = "fixed") +
    labs(x = "t", y = "Value",
         title = "Diff  ")



  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Example time grid
  t <- seq(0, 10, length.out = 100)

  # Example functions for b_{ij}(t)
  b10 <- sin(t)
  b20 <- cos(t)
  b30 <- t^0.5


  b11 <- 0.2 * t
  b21 <- -0.1 * t
  b31 <- 0.05 * t

  # Compute probabilities
  exp1 <- exp(b10)
  exp2 <- exp(b20)
  exp3 <- exp(b30)
  p1 <- exp1 / (exp1 + exp2 + exp3)

  exp1_shift <- exp(b10 + b11)
  exp2_shift <- exp(b20 + b21)
  exp3_shift <- exp(b30 + b31)
  p2 <- exp1_shift / (exp1_shift + exp2_shift + exp3_shift)

  # Difference
  delta <- p2 - p1

  # Data for plotting
  df <- data.frame(t = t, Difference = delta, p1 = p1, p2 = p2)

  # Plot
  ggplot(df, aes(x = t, y = Difference)) +
    geom_line(color = "blue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Difference in Probability for Class 1 Over Time",
         y = expression(Delta(t) == p[2](t) - p[1](t)),
         x = "Time") +
    theme_minimal()


  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)  # optional for multi-panel layout

  # Time grid
  t <- seq(0, 1, length.out = 100)

  # Example functions for b_ij(t)
  b10 <-  2 * sin(pi * t) - 1
  b20 <-  2 * sin(pi * t) - 1
  b30 <- 2 - 2 * t

  b11 <- t
  b21 <- -t
  b31 <- t

  # Compute baseline probabilities
  exp1 <- exp(b10)
  exp2 <- exp(b20)
  exp3 <- exp(b30)


  # Compute slope-adjusted probabilities
  exp1_shift <- exp(b10 + b11)
  exp2_shift <- exp(b20 + b21)
  exp3_shift <- exp(b30 + b31)

  p1_0 <- exp1 / (exp1 + exp2 + exp3)
  p1_1 <- exp1_shift / (exp1_shift + exp2_shift + exp3_shift)

  p2_0 <- exp2 / (exp1 + exp2 + exp3)
  p2_1 <- exp2_shift / (exp1_shift + exp2_shift + exp3_shift)

  p3_0 <- exp3 / (exp1 + exp2 + exp3)
  p3_1 <- exp3_shift / (exp1_shift + exp2_shift + exp3_shift)

  # Compute difference
  delta1 <- p1_0 - p1_1
  delta2 <- p2_0 - p2_1
  delta3 <- p3_0 - p3_1


  # Response 1
  # Combine into a tidy data frame
  df <- data.frame(t = t, Z0 = p1_0, Z1 = p1_1, Difference = delta1) %>%
    pivot_longer(cols = c("Z0", "Z1"), names_to = "Curve", values_to = "Probability")

  # Plot 1: Original curves
  p_curves <- ggplot(df, aes(x = t, y = Probability, color = Curve)) +
    geom_line(size = 1) +
    labs(title = "exp (Generating curves) - Response 1",
         y = "Alpha", x = "Time")
  # Plot 2: Difference
  p_diff <- ggplot(df %>% filter(Curve == "Z0") %>%
                     mutate(Difference = delta1),
                   aes(x = t, y = Difference)) +
    geom_line(color = "blue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Difference in Class 1 Probability (Slope - Baseline)",
         y = "y", x = "Time")

  # Combine plots (stacked vertically)
  r1 <- p_curves / p_diff

  # Response 2
  # Combine into a tidy data frame
  df <- data.frame(t = t, Z0 = p2_0, Z1 = p2_1, Difference = delta2) %>%
    pivot_longer(cols = c("Z0", "Z1"), names_to = "Curve", values_to = "Probability")

  # Plot 1: Original curves
  p_curves <- ggplot(df, aes(x = t, y = Probability, color = Curve)) +
    geom_line(size = 1) +
    labs(title = "exp (Generating curves) - Response 1",
         y = "Alpha", x = "Time")
  # Plot 2: Difference
  p_diff <- ggplot(df %>% filter(Curve == "Z0") %>%
                     mutate(Difference = delta2),
                   aes(x = t, y = Difference)) +
    geom_line(color = "blue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Difference in Class 1 Probability (Slope - Baseline)",
         y = "y", x = "Time")

  # Combine plots (stacked vertically)
  r2 <- p_curves / p_diff

  # Response 2
  # Combine into a tidy data frame
  df <- data.frame(t = t, Z0 = p3_0, Z1 = p3_1, Difference = delta3) %>%
    pivot_longer(cols = c("Z0", "Z1"), names_to = "Curve", values_to = "Probability")

  # Plot 1: Original curves
  p_curves <- ggplot(df, aes(x = t, y = Probability, color = Curve)) +
    geom_line(size = 1) +
    labs(title = "exp (Generating curves) - Response 1",
         y = "Alpha", x = "Time")
  # Plot 2: Difference
  p_diff <- ggplot(df %>% filter(Curve == "Z0") %>%
                     mutate(Difference = delta3),
                   aes(x = t, y = Difference)) +
    geom_line(color = "blue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Difference in Class 1 Probability (Slope - Baseline)",
         y = "y", x = "Time")

  # Combine plots (stacked vertically)
  r3 <- p_curves / p_diff

  r1 + r2 + r3 + plot_layout(ncol = 3)
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


