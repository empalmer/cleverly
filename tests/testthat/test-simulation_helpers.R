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


test_that("No Z - Functional simulation", {
  testthat::skip("Debugging, not test")

  # Generate simulation data
  set.seed(123)
  sim <- sim_noZ()
  time <- sim$time
  id <- sim$individual
  time_id <- sim$Capture.Number
  Yij0 <- sim$total_n
  Y <- dplyr::select(sim, -c(
                        "individual",
                        "time",
                        "total_n",
                        "Capture.Number"))
  Y <- as.matrix(Y)
  K <- ncol(Y)

  # Run algorithm
  iter <- 25
  res <- cleverly(Y = Y,
                  subject_ids = id,
                  time = time,
                  gammas = 1000, # controls smoothness
                  tau = 1/2,
                  theta = 300,
                  psi = 8000, # controls clustering
                  C = 100,
                  max_admm_iter = iter,
                  max_outer_iter = 1)


  # Diagnostics:
  phis <- res$result$phis_list[[1]]
  v <- res$result$v
  v
  v_mat <- matrix(v, nrow = 6)

  # See the final fit
  visualize_final_fit(Y = Y,
                      res = res)

  # View how curves change over iterations.
  betas <- res$result$admm_beta_list[[1]]
  beta_path(betas,
            K = K,
            time = time)

  # Examine differences:
  # Differences in the ADMM iterations
  diffs <- res$result$admm_diffs[[1]]
  plot_differneces(diffs)
  diffs

  # Examine clusters:
  A <- res$result$clusters

})





test_that("Bspline sim no Z",{
  testthat::skip("Debugging, not test")
  set.seed(123)
  n <- 100
  time_list <- sim_timepoints(n = n)
  time <- time_list$X$time
  mi_vec <- time_list$mi_vec
  is <- rep(1:n, mi_vec)
  M <- sum(mi_vec)
  K <- 4
  Z <- matrix(rep(1, M), nrow = M)
  B <- get_B(time, order = 3, nknots = 3)
  P <- 6







  # 2 clusters
  betaC1 <- matrix(c(rep(c(1, 2, 3, 4, 5, 6), 2), ncol = 1))#l = 2
  betaC2 <- matrix(.5*rev(betaC1), ncol = 1)
  beta <- rbind(betaC1, betaC2)

  Y <- sim_Y(beta = beta,
             Z = Z,
             B = B,
             K = K,
             mi_vec = mi_vec)
  Y_ra <- Y/rowSums(Y)

  # Plot data:
  data.frame(time = time, Y_ra) %>%
    tidyr::pivot_longer(-time) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~name)

  # Truth:
  y_hat_truth <- estimate_y(beta = beta,
                       B = B,
                       Z = Z,
                       K = K)
  data.frame(time = time, y_hat_truth) %>%
    tidyr::pivot_longer(-time) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~name)


  # Run code
  iter <- 25
  res <- cleverly(Y = Y,
                  subject_ids = is,
                  time = time,
                  gammas = .1,
                  psi = .01,
                  max_admm_iter = iter,
                  max_outer_iter = 1)

  visualize_curve(beta_new, B, Z, K = 4, time)

  betas <- res$result$admm_beta_list[[1]]
  plot <- beta_path(betas, K)


})



test_that("Simulation With Z", {

  skip("Skip")
  # Generate simulation data
  set.seed(127)
  sim <- sim_Z_longitudinal(n = 20,
                            range_start = 5000,
                            range_end = 10000,
                            nknots = 3,
                            K = 12,
                            order = 3,
                            user_var = 1000,
                            cor_str = "CON",
                            rho = 0.8,
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
                                 color = factor(Z))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::facet_wrap(~name) +
    ggplot2::labs(title = "Simulated Data",
                  color = "EV",
                  y = "Count",
                  x = "Time")



  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))
  Z <- sim$Z


  tau = .005
  theta = 300
  psi = 800
  gammas = c(10, 10)
  #start <- Sys.time()
  #Rprof("test.out", interval = .02)
  profvis::profvis({
  res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  lp = 0,
                  time = time,
                  # Hyperparameters
                  gammas = gammas,
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
                  cor_str = "CON")
  })
  #end <- Sys.time()
  #Rprof(NULL)
  #summaryRprof("test.out")$by.self[1:10,1:2]
  #(duration <- end - start)
  res$clusters

  # Diagnostic plots:
  plot_clusters(res = res,
                K = 12,
                tau = 0,
                psi = 0,
                theta = 0,
                gammas = 0,
                max_admm_iter = 0,
                max_outer_iter = 0)

  plot_cluster_path(res, psi, tau, theta,gammas, max_admm_iter, max_outer_iter, duration)
  plot_initial_fit(res, K = 12, gammas = gammas)

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


})



test_that("psi BIC", {

  skip("Skip")

  set.seed(123)
  sim <- sim_Z_longitudinal()
  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))
  Z <- sim$Z


  start <- Sys.time()
  psis <- seq(3000, 3300, by = 75)
  gammas <- c(5000, 10000)
  hyps <- expand.grid(psi = psis,
                      gammas = gammas)
  res_list <- purrr::pmap(hyps, ~cleverly(Y = Y,
                                         Z = Z,
                                         subject_ids = individual,
                                         lp = 0,
                                         time = time,
                                         gammas = c(.1*..2, ..2), # controls smoothness
                                         tau = .1, # Controls cuttoff shrinkage
                                         theta = 3000, # for lambda, but also for d
                                         psi = ..1, # controls clustering
                                         C = 100,
                                         max_admm_iter = 100,
                                         max_outer_iter = 10,
                                         max_2_iter = 100,
                                         epsilon_r = .001,
                                         epsilon_d = .05,
                                         epsilon_b = .001,
                                         epsilon_2 = .001))
  end <- Sys.time()
  (duration <- end - start)


  psis <- c(200, 400, 500)
  taus <- c(.001, .01, .05)
  thetas <- c(150, 200)
  max_admm_iter = 100
  max_outer_iter = 60
  hyps <- expand.grid(psi = psis,
                      taus = taus,
                      theta = thetas)
  start <- Sys.time()
  future::plan(multisession, workers = 8)
  res_list <- furrr::future_pmap(hyps, ~cleverly(Y = Y,
                                                 Z = Z,
                                                 subject_ids = individual,
                                                 lp = 0,
                                                 time = time,
                                                 gammas = c(5,5),
                                                 tau = ..2,
                                                 theta = ..3,
                                                 psi = ..1,
                                                 C = 100,
                                                 max_admm_iter = max_admm_iter,
                                                 max_outer_iter = max_outer_iter,
                                                 epsilon_r = 1,
                                                 epsilon_d = 1,
                                                 epsilon_b = 1))
  end <- Sys.time()
  (duration <- end - start)


  hyps[best,]
  psi <- hyps[best,1]
  tau <- hyps[best,2]
  theta <- hyps[best, 3]


  purrr::map_dbl(res_list, ~.x$BIC)
  best <- which.min(purrr::map_dbl(res_list, ~.x$BIC))

  res <- res_list[[best]]


  purrr::map_dbl(res_list, ~.x$BIC)
  (psi <- psis[best])

  y_hat <- res$y_hat
ß
  res <- res_list[[5]]

  # D convergence
  purrr::imap_dfr(res$d_list,
                  ~data.frame(cluster = unlist(.x),
                              run = .y)) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(x = row, y = cluster, color = factor(run))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 3, linetype = 2) +
    ggplot2::labs(x = "Overall Iteration",
                  y = "Number of Clusters",
                  title = "d",
                  subtitle = paste("psi:",psi,
                                   ",tau:",round(tau,2),
                                   "theta:",theta,
                                   ",admm_iter:",max_admm_iter,
                                   "outer_iter:",max_outer_iter,
                                   "time:", duration),
                  color = "Outer iteration")
  # r convergence
  purrr::imap_dfr(res$r_list,
                  ~data.frame(cluster = unlist(.x),
                              run = .y)) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(x = row, y = cluster, color = factor(run))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 3, linetype = 2) +
    ggplot2::labs(x = "Overall Iteration",
                  y = "Number of Clusters",
                  title = "r = beta_k - beta_k' - v_kappa",
                  subtitle = paste("\u03A8: ", psi,
                                   ",\u03C4:",round(tau,2),
                                   ",\u03B8:",theta,
                                   "gamma: ", gamma,
                                   ",admm",max_admm_iter,
                                   "outer",max_outer_iter,
                                   "time:", duration),
                  color = "Outer iteration")

  # Cluster progress:
  cluster_track <- res$cluster_list
  purrr::imap_dfr(cluster_track,
                  ~data.frame(cluster = purrr::map_dbl(.x, ~.x$no),
                              run = .y)) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(x = row, y = cluster, color = factor(run))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 3, linetype = 2) +
    ggplot2::labs(x = "Overall Iteration",
                  y = "Number of Clusters",
                  title = "Cluster Progress",
                  subtitle = paste("psi:",psi,
                                   ",tau:",round(tau,2),
                                   "theta:",theta,
                                   ",admm_iter:",max_admm_iter,
                                   "outer_iter:",max_outer_iter,
                                   "time:", duration),
                  color = "Outer iteration")

  # Get cluster membership
  cluster_df <- data.frame(
    K = factor(1:12, levels = 1:12),
    cluster = factor(res$clusters$membership))
  knitr::kable(table(cluster_df$cluster, rep(1:3, each = 4)))

  # plot clusters
  y_hat %>%
    dplyr::mutate(response = factor(response, levels = 1:12)) %>%
    dplyr::left_join(cluster_df, by = c("response" = "K")) %>%
    dplyr::mutate(clusterZ = ifelse(Z == 1, "EV", cluster)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time)) +
    ggplot2::geom_point(ggplot2::aes(y = y,
                                     color = factor(Z),
                                     shape = factor(Z)),
                        size = .6, alpha = .8) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_line(ggplot2::aes(y = yhat,
                                    color = clusterZ,
                                    group = factor(Z)),
                       linewidth = 1) +
    ggplot2::facet_wrap(~response) +
    ggplot2::scale_color_manual(
      # values = c(RColorBrewer::brewer.pal(n = length(unique(cluster_df$cluster)),
      #                       name = "Set1"),
      #            "grey50"),
      #values = c(rcartocolor::carto_pal(length(unique(cluster_df$cluster)),
      #                                  "Safe"),
      #           "grey50"),
      values = c(viridis::viridis(length(unique(cluster_df$cluster))), "grey50"),
      name = "Cluster"
    ) +
    ggplot2::labs(subtitle   = paste0("psi:",psi,
                                      ", tau:",round(tau,2),
                                      ", theta:",theta,
                                      ", admm_iter:",max_admm_iter,
                                      ", outer_iter:",max_outer_iter))

})





test_that("cleverly best psi", {

  skip("Skip")
  # Generate simulation data
  set.seed(127)
  sim <- sim_Z_longitudinal(n = 20,
                            range_start = 4000,
                            range_end = 8000,
                            nknots = 3,
                            K = 12,
                            order = 3,
                            user_var = 600,
                            cor_str = "CON",
                            rho = 0.4,
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
                                 color = factor(Z))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::facet_wrap(~name, scales = "free_y") +
    ggplot2::labs(title = "Simulated Data",
                  color = "EV",
                  y = "Count",
                  x = "Time")



  Y <- dplyr::select(sim, -c(
    "total_n",
    "Capture.Number",
    "Z"))
  Z <- sim$Z

  res_psi <- cleverly_bestpsi(psi_min = 800,
                          psi_max = 1000,
                          npsi = 2,
                          parralel = FALSE,
                          Y = Y,
                          Z = Z,
                          lp = 0,
                          time = time,
                          # Hyperparameters
                          gammas = c(1, 1),
                          tau = 0.005,
                          theta = 300,
                          C = 100,
                          # Iterations max
                          max_admm_iter = 50,
                          max_outer_iter = 2,
                          max_2_iter = 50,
                          # Convergence criteria
                          epsilon_r = .001,
                          epsilon_d = .05,
                          epsilon_b = .01,
                          epsilon_2 = .001,
                          cor_str = "IND")

  res_psi$clusters

  res <- res_psi

  # Diagnostic plots:
  plot_clusters(res = res,
                K = 12,
                tau = tau,
                psi = psi,
                theta = theta,
                gammas = gammas,
                max_admm_iter = max_admm_iter,
                max_outer_iter = max_outer_iter)

  plot_cluster_path(res, psi, tau, theta,gammas, max_admm_iter, max_outer_iter, duration)
  plot_initial_fit(res, K = 12, gammas = gammas)

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


})


test_that("cleverly best psi", {
  set.seed(127)

  sim <- one_sim()

  full_sim <- my_method_sim(nsim = 2)
  full_sim2 <- my_method_sim(nsim = 2)
  purrr::map_dfr(full_sim, "cluster_result")
  purrr::map(full_sim, "all_clusters_psi")


  full_sim3 <- my_method_sim(nsim = 1, parralel = T)

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


test_that("Simulation filter Z", {

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
                            al = 0.4,
                            slope_base = "cluster_base_alldiff_slope")

  sim_base_only <- sim %>%
    dplyr::filter(Z == 0)

  # Visualize simulated data
  sim_base_only %>%
    tidyr::pivot_longer(-c(individual,
                           time,
                           Capture.Number,
                           total_n, Z)) %>%
    dplyr::mutate(name = factor(name,
                                levels = paste0("Taxa.", 1:12))) %>%
    ggplot2::ggplot(ggplot2::aes(x = time,
                                 y = value)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::facet_wrap(~name) +
    ggplot2::labs(title = "Simulated Data",
                  color = "EV",
                  y = "Count",
                  x = "Time")



  Y <- dplyr::select(sim_base_only, -c(
    "total_n",
    "Capture.Number",
    "Z"))

  tau = .01
  theta = 300
  psi = 500
  gammas = c(1,1)
  start <- Sys.time()
  #Rprof("test.out", interval = .02)
  #profvis::profvis({
  res <- cleverly(Y = Y,
                  subject_ids = individual,
                  lp = 0,
                  time = time,
                  # Hyperparameters
                  gammas = c(1),
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
  res$clusters$no


  res_psi <- cleverly_bestpsi(psi_min = 100,
                              psi_max = 2000,
                              npsi = 2,
                              parralel = F,
                              Y = Y,
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


  # Diagnostic plots:
  plot_clusters(res = res,
                K = 12,
                tau = 0,
                psi = 0,
                theta = 0,
                gammas = 0,
                max_admm_iter = 0,
                max_outer_iter = 0)

  plot_cluster_path(res, psi, tau, theta,gammas, max_admm_iter, max_outer_iter, duration)
  plot_initial_fit(res, K = 12, gammas = gammas)

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


})
