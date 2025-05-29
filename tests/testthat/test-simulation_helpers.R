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
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/cond_large_rho_var/sim_data_2.rds")
  sim <- read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/cond_small_var_Y0/sim_data_2.rds")
  sim <- read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/cond_1mil_similar_ranges/sim_data_2.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/CONd_may22_1mil/sim_data_74.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_binary_baseline/sim_data_27.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_binary_slope/sim_data_27.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_cont_baseline/sim_data_1.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_cont_slope/sim_data_27.rds")

  # Visualize simulated data
  Z_type = "continuous"
  Z_type = "binary"
  plot_sim_data(sim, Z_type = Z_type)

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
                  lp = 0,
                  cor_str = "IND",
                  # Hyperparameters
                  gammas = c(1,1),
                  theta = 500,
                  parralel = F,
                  psi_min = 50,
                  psi_max = 228.57,
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
                Z_type = Z_type,
                y_type = "y_hat_baseline")

  plot_clusters(res,
                response_names = 1:12,
                y_type = "slope",
                Z_type = Z_type,
                scales = "fixed")


  plot_clusters(res,
                Z = rep(Z, 12),
                response_names = 1:12,
                Z_type = "binary",
                y_type = "y_hat_counts_group")


  plot_BIC(res, BIC_type = "BIC", psis = seq(600, 1400, length.out = 2))
  plot_BIC(res, BIC_type = "BIC_ra_group", psis = seq(600, 1400, length.out = 2))

  res$possible_clusters
  res$psi
  res$phi
  res$rho

  res$BIC

  res$clusters



  plot_initial_fit(res, K = 12)

  # Get cluster membership
  cluster_df <- data.frame(
    K = factor(1:12, levels = 1:12),
    cluster = factor(res$clusters$membership))
  knitr::kable(table(cluster_df$cluster, rep(1:3, each = 4)))


  get_corR("CON", mi = 3, K = 2, rho = .4)


})



test_that("CLR", {
  skip("Skip - used as test file for cleverly")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/cond_large_rho_var/sim_data_2.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/cond_small_var_Y0/sim_data_2.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/cond_1mil_similar_ranges/sim_data_2.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/CONd_may22_1mil/sim_data_74.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_binary_baseline/sim_data_27.rds")
  sim <- readr::read_rds("~/Desktop/Research/buffalo-sciris/novus_results/sim_data/Z_cont_baseline/sim_data_27.rds")


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
                         M = nrow(Y)) %>%
    get_cluster_diagnostics(true_cluster)
  clr_res$cluster_diagnostics



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


test_that("real data", {

  skip("Skip")

#
#   library(phyloseq)
#   library(tidyverse)

  ps <- readRDS("~/Desktop/Research/buffalo-sciris/buffalo/Buffalo_Data/phyloseq_to_publish.rds")



  total_gen <- otu_table(ps) %>%
    base::colMeans() %>%
    sort()
  length(names(total_gen))
  selected_genera <- names(total_gen[total_gen > 89.75])
  length(selected_genera)



  ps_filtered <- subset_taxa(ps, Genus %in% selected_genera)
  ps_filtered

  metadata <- sample_data(ps_filtered)
  names(metadata)


  # Extract informatino needed for clevelry

  time <- as.numeric(metadata$Capture)

  time <- as.numeric(metadata$Age_At_Capture)
  # Check that this code is working...
  ids <- as.numeric(metadata$AnimalID)
  Z <- as.numeric(metadata$High_Nutrition) - 1

  Y <- otu_table(ps_filtered, taxa_are_rows = FALSE) %>%
    data.frame()

  library(cleverly)
  res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = ids,
                  lp = 0,
                  time = time,
                  gammas = c(1, 1), # controls smoothness
                  tau = 0.005, # Controls cuttoff for highest shrinkage
                  theta = 300, # for lambda, but also for d
                  psi = 800, # controls clustering
                  C = 100,
                  max_admm_iter = 100,
                  max_outer_iter = 10,
                  max_2_iter = 100,
                  cor_str = "IND")

res$clusters
plot_clusters(res = res,
                K = 32,
              EV = "grey60")




cluster_df <- data.frame(
  K = factor(1:K, levels = 1:K, labels = names(Y)),
  cluster = factor(res$clusters$membership)) %>%
  arrange(cluster, K) %>%
  mutate(K = factor(K, levels = K, ))

K = 32
ordering <- cluster_df$K

# plot clusters
res$y_hat %>%
  dplyr::mutate(response = factor(response, levels = 1:K, labels = names(Y))) %>%
  dplyr::left_join(cluster_df, by = c("response" = "K")) %>%
  dplyr::mutate(clusterZ = ifelse(Z == 1, "EV", cluster),
                response = factor(response, levels = ordering, ordered = T)) %>%
  ggplot2::ggplot(ggplot2::aes(x = time)) +
  ggplot2::geom_point(ggplot2::aes(y = y,
                                   color = factor(Z),
                                   shape = factor(Z)),
                      size = 1, alpha = .8) +
  ggplot2::labs(color = "EV",
                shape = "EV") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_line(ggplot2::aes(y = yhat,
                                  color = clusterZ,
                                  group = factor(Z)),
                     linewidth = 1) +
  ggplot2::facet_wrap(~response, scales = "free_y", nrow = 4) +
  ggplot2::scale_color_manual(
    values = c(sample(viridis::viridis(length(unique(cluster_df$cluster)))), "grey"),
    name = "Cluster"
  )

# Start 1053
res_list <- list()
psis <- seq(500, 1500, length.out = 5)
for (p in 1:length(psis)) {
  psi <- psis[p]
  print(psi)
  res_list[[p]] <- cleverly(Y = Y,
                            Z = Z,
                            subject_ids = ids,
                            time = time,
                            lp = 0,
                            cor_str = "IND",
                            max_outer_iter = 10,
                            max_admm_iter = 100,
                            max_2_iter = 200,
                            gammas = c(1, 1), # controls smoothness
                            tau = 0.005, # Controls cuttoff for highest shrinkage
                            theta = 300, # for lambda, but also for d
                            psi = psi) # controls clustering)
}

best <- which.min(purrr::map_dbl(res_list, ~.x$BIC))
res <- res_list[[best]]

})
