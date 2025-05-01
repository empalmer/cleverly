test_that("timepoints work", {
  time <- sim_timepoints(n = 5)
  expect_type(time, "list")
})

test_that("Z works", {
  time <- sim_timepoints(n = 5)
  mi_vec <- time$mi_vec
  Z_sim <- sim_Z(mi_vec)
})

test_that("Test sim Y works", {
  set.seed(124)
  time_list <- sim_timepoints(n = 5)
  mi_vec <- time_list$mi_vec
  i_index <- c(0, cumsum(mi_vec))
  time <- time_list$X$time

  K <- 4
  Z <- sim_Z(mi_vec)

  B <- get_B(time, order = 3, nknots = 3)

  # 2 clusters
  betaC1 <- matrix(c(rep(c(1, 1, 1, 1, 1, 1), 2),     #l = 0
                     rep(c(-2, -2, -2, -2, -2, -2),2),  #l = 1
                     rep(c(1, 2, 3, 4, 5, 6), 2)), ncol = 3)  #l = 2
  betaC2 <- .5*betaC1

  beta <- rbind(betaC1, betaC2)

  Y_ij <- sim_Y_ij(i = 1,
                   j = 1,
                   beta,
                   Z,
                   B,
                   Y_ij0 = 100,
                   K = K,
                   mi_vec = mi_vec,
                   i_index = i_index)

  Y_i <- sim_Yi(i = 1,
                beta = beta,
                Z = Z,
                B = B,
                K = K,
                mi_vec = mi_vec,
                i_index = i_index)

  Y <- sim_Y(beta = beta,
             Z = Z,
             B = B,
             K = K,
             mi_vec = mi_vec,
             i_index = i_index)
})



test_that("Simulation With Z", {
  skip("Skip")
  # Generate simulation data
  set.seed(127)
  sim <- sim_Z_longitudinal(n = 20,
                            range_start = 5000,
                            range_end = 20000,
                            nknots = 3,
                            K = 12,
                            order = 3,
                            user_var = 10000,
                            cor_str = "CON-d",
                            rho = 0.6,
                            prob1 = .5,
                            slope_base = "cluster_base_alldiff_slope")

  # Visualize simulated data
  sim %>%
    tidyr::pivot_longer(-c(individual,
                           time,
                           Capture.Number,
                           total_n, Z)) %>%
    dplyr::mutate(name = factor(name,
                                levels = paste0("Taxa.", 1:12))) %>%
    ggplot2::ggplot(ggplot2::aes(x = time,
                                 y = value,
                                 color = factor(Z),
                                 shape = factor(Z))) +
    ggplot2::geom_jitter(size = 1) +
    ggplot2::facet_wrap(~name) +
    ggplot2::labs(title = "Simulated Data",
                  color = "EV",
                  shape = "EV",
                  y = "Count",
                  x = "Time")



  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))
  Z <- sim$Z

  psi = 800
  #start <- Sys.time()
  #Rprof("test.out", interval = .02)
  #profvis::profvis({
  res <- cleverly_onepsi(Y = Y,
                         Z = Z,
                         subject_ids = individual,
                         time = time,
                         lp = 0,
                         cor_str = "AR1-d",
                         # Hyperparameters
                         gammas = c(1,1),
                         psi = psi,
                         # Iterations max
                         max_admm_iter = 10,
                         max_outer_iter = 1,
                         max_2_iter = 10,
                         )
  #})
  # end <- Sys.time()
  # Rprof(NULL)
  # summaryRprof("test.out")$by.self[1:15,1:2]
  #(duration <- end - start)
  res$clusters

  # Diagnostic plots:
  plot_clusters(res = res,
                response_names = LETTERS[1:12])


  plot_initial_fit(res, K = 12)

  # Get cluster membership
  cluster_df <- data.frame(
    K = factor(1:12, levels = 1:12),
    cluster = factor(res$clusters$membership))
  knitr::kable(table(cluster_df$cluster, rep(1:3, each = 4)))

  res_psi <- cleverly(Y = Y,
                      Z = Z,
                      subject_ids = individual,
                      time = time,
                      lp = 0,
                      cor_str = "AR1-d",
                      # Hyperparameters
                      gammas = c(1,1),
                      npsi = 2,
                      # Iterations max
                      max_admm_iter = 10,
                      max_outer_iter = 1,
                      max_2_iter = 10,
  )


})




test_that("Simulation no Z", {

  skip("Skip")
  # Generate simulation data
  set.seed(127)
  sim <- sim_Z_longitudinal(n = 20,
                            range_start = 5000,
                            range_end = 10000,
                            nknots = 3,
                            K = 12,
                            order = 3,
                            user_var = 500,
                            cor_str = "IND",
                            rho = 0.4,
                            slope_base = "cluster_base_alldiff_slope")

  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))


  tau = .01
  theta = 300
  psi = 500
  start <- Sys.time()
  #Rprof("test.out", interval = .02)
  #profvis::profvis({
  res_noz <- cleverly(Y = Y,
                  subject_ids = individual,
                  lp = 0,
                  time = time,
                  # Hyperparameters
                  gammas = c(10),
                  tau = tau,
                  theta = theta,
                  psi = psi,
                  C = 100,
                  # Iterations max
                  max_admm_iter = 100,
                  max_outer_iter = 10,
                  max_2_iter = 100,
                  # Convergence criteria
                  epsilon_r = .001,
                  epsilon_d = .05,
                  epsilon_b = .01,
                  epsilon_2 = .001,
                  cor_str = "IND")
  # })
  end <- Sys.time()
  #Rprof(NULL)
  #summaryRprof("test.out")$by.self[1:10,1:2]
  (duration <- end - start)
  res_noz$clusters$no

  cluster_df <- data.frame(
    K = factor(1:K, levels = 1:K),
    cluster = factor(res$clusters$membership))

  y_hat <- res_noz$y_hat
  # plot clusters
  y_hat %>%
    dplyr::mutate(response = factor(response, levels = 1:K)) %>%
    dplyr::left_join(cluster_df, by = c("response" = "K")) %>%
    ggplot2::ggplot(ggplot2::aes(x = time)) +
    ggplot2::geom_point(ggplot2::aes(y = y),
                        size = .6, alpha = .8) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_line(ggplot2::aes(y = yhat,
                                    color = cluster),
                       linewidth = 1) +
    ggplot2::facet_wrap(~response)


  # Diagnostic plots:
  plot_clusters(res = res_noz,
                K = 12,
                tau = tau,
                psi = psi,
                theta = 0,
                gammas = 1,
                max_admm_iter = 100,
                max_outer_iter = 0)

  plot_cluster_path(res,0, 0, 0,0,0, 0, 0)
  plot_initial_fit(res_noz, K = 12, gammas = 1)

  plot_alg2_convergence(res)
  plot_d_convergence(res)
  plot_r_convergence(res)

  # Get cluster membership
  cluster_df <- data.frame(
    K = factor(1:12, levels = 1:12),
    cluster = factor(res$clusters$membership))
  knitr::kable(table(cluster_df$cluster, rep(1:3, each = 4)))

  # number of iterations
  unlist(res$ts) #ADMM
  unlist(res$rs) #alg2
  # convergence
  unlist(res$d_list)
  unlist(res$r_list)
  unlist(res$alg_2_beta_diff)
  unlist(res$alg1_diff)


  y_hat_init <- res_noz$y_hat_init
  # plot clusters

  y_hat <- estimate_y(beta = beta,
                      B = B,
                      Z = Z,
                      K = K,
                      Y = Y,
                      time = time)
  y_hat %>%
    dplyr::mutate(response = factor(response, levels = 1:K)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time)) +
    ggplot2::geom_point(ggplot2::aes(y = y),
                        size = .6, alpha = .8) +
    ggplot2::geom_line(ggplot2::aes(y = yhat),
                       linewidth = 1, color = "blue") +
    ggplot2::facet_wrap(~response)


  res_psi <- cleverly_bestpsi(psi_min = 100,
                              psi_max = 2000,
                              npsi = 2,
                              parralel = T,
                              Y = Y,
                              Z = matrix(1, nrow = nrow(Y), ncol = 1),
                              lp = 0,
                              time = time,
                              # Hyperparameters
                              gammas = c(10),
                              tau = 0.01,
                              theta = 300,
                              C = 100,
                              # Iterations max
                              max_admm_iter = 200,
                              max_outer_iter = 10,
                              max_2_iter = 100,
                              # Convergence criteria
                              epsilon_r = .001,
                              epsilon_d = .05,
                              epsilon_b = .01,
                              epsilon_2 = .001,
                              cor_str = "IND")
  res_psi <- cleverly_bestpsi(psi_min = 100,
                              psi_max = 2000,
                              npsi = 2,
                              parralel = F,
                              Y = Y,
                              Z = matrix(1, nrow = nrow(Y), ncol = 1),
                              lp = 0,
                              time = time,
                              # Hyperparameters
                              gammas = c(10),
                              tau = 0.01,
                              theta = 300,
                              C = 100,
                              # Iterations max
                              max_admm_iter = 200,
                              max_outer_iter = 10,
                              max_2_iter = 100,
                              # Convergence criteria
                              epsilon_r = .001,
                              epsilon_d = .05,
                              epsilon_b = .01,
                              epsilon_2 = .001,
                              cor_str = "IND")

})



test_that("real data", {

  skip("Skip")


  library(phyloseq)
  library(tidyverse)

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

}
