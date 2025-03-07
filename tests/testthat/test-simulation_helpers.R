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
                   mi_vec = mi_vec)

  Y_i <- sim_Yi(i = 1,
                beta = beta,
                Z = Z,
                B = B,
                K = K,
                mi_vec = mi_vec)

  Y <- sim_Y(beta = beta,
             Z = Z,
             B = B,
             K = K,
             mi_vec = mi_vec)
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
  set.seed(123)
  sim <- sim_Z_longitudinal()
  sim %>%
    tidyr::pivot_longer(-c(individual,
                           time,
                           Capture.Number,
                           total_n, Z)) %>%
    dplyr::mutate(name = factor(name,
                         levels = paste0("Taxa.", 1:12) )) %>%
    ggplot2::ggplot(ggplot2::aes(x = time,
                                 y = value,
                                 color = factor(Z))) +
    ggplot2::geom_point(size = 1) +
    ggplot2::facet_wrap(~name)

  Y <- dplyr::select(sim, -c(
                        "total_n",
                        "Capture.Number",
                        "Z"))
  Z <- sim$Z


  psi <- 200
  tau <- 1/2
  theta <- 300
  max_admm_iter = 20
  max_outer_iter = 2
  start <- Sys.time()
  Rprof("test.out")
  res <- cleverly(Y = Y,
                  Z = Z,
                  subject_ids = individual,
                  lp = 0,
                  time = time,
                  gammas = c(100,100), # controls smoothness
                  tau = tau, # Controls cuttoff for highest shrinkage
                  theta = theta, # for lambda, but also for d
                  psi = psi, # controls clustering
                  C = 100,
                  max_admm_iter = max_admm_iter,
                  max_outer_iter = max_outer_iter,
                  epsilon_r = 1,
                  epsilon_d = 1,
                  epsilon_b = 1)
  end <- Sys.time()
  Rprof(NULL)
  summaryRprof("test.out")$by.self[1:10,]
  duration <- end - start

  y_hat <- res$y_hat

  # Examine convergence:
  # Algorithm 1 iterations
  (s <- length(res$admm_beta_list))
  # Algorithm 3 admm iterations
  (t <- purrr::map_dbl(res$admm_beta_list, length))

  # admm parts:
  purrr::map(res$r_list, unlist) # beta - beta - v
  purrr::map(res$d_list, unlist) # sum of vs
  # backfitting parts:
  unlist(res$loop_list_diff)

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
                  subtitle = paste("psi:",psi,
                                   ",tau:",round(tau,2),
                                   "theta:",theta,
                                   ",admm_iter:",max_admm_iter,
                                   "outer_iter:",max_outer_iter,
                                   "time:", duration),
                  color = "Outer iteration")

  # Cluster progress:
  cluster_track <- res$cluster_list
  purrr::imap_dfr(cluster_track,
                           ~data.frame(cluster = purrr::map_dbl(.x, ~.x$no),
                                       run = .y)) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(x = row, y = cluster, color = factor(run))) +
    ggplot2::geom_point() +
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
