# Plots -------------------------------------------------------------------

#' Plot the final fit
#'
#' Plot the final fit of the model (relative abundance scale) with the data points (relative abundance scale)
#'
#' @param Y Input matrix Matrix of counts Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param res cleverly model output object
#'
#' @returns ggplot object
#' @export
visualize_final_fit <- function(Y, res){
  K <- ncol(Y)
  beta <- res$result$beta

  # Rename yhats
  y_hat <- res$result$y_hat
  colnames(y_hat) <- paste0("Response.", 1:K)

  # Convert Y s to RA
  y_ra <- Y/rowSums(Y)
  colnames(y_ra) <- paste0("Response.", 1:K)

  y_ra <- data.frame(y_ra) %>%
    dplyr::mutate(time = time) %>%
    tidyr::pivot_longer(-time,
                        names_to = "response",
                        values_to = "Y")
  y_hat_ra <- data.frame(y_hat) %>%
    dplyr::mutate(time = time) %>%
    tidyr::pivot_longer(-time,
                        names_to = "response",
                        values_to = "y_hat")

  # Get cluster membership
  cluster_df <- data.frame(response = paste0("Response.", 1:K),
                           cluster = factor(sim_res$result$clusters$membership))

  # Y-hats
  Y_ra_df <- y_ra %>%
    dplyr::mutate(y_hat = y_hat_ra$y_hat) %>%
    dplyr::left_join(cluster_df,
                     by = c("response" = "response")) %>%
    dplyr::mutate(taxa = factor(response,
                                levels = paste0("Response.", 1:K)))

  # Make plot
  plot <- Y_ra_df %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = time, y = Y),
                        size = 1, color = "black") +
    ggplot2::geom_line(ggplot2::aes(x = time,
                                    y = y_hat,
                                    color = cluster),
                       linewidth = 2) +
    ggplot2::facet_wrap(~response) +
    ggplot2::labs(title = paste0("Model fit (penalized"),
                  x = "Time",
                  y = "Relative abundances",
                  color = "Cluster membership")

  return(plot)

}

#' Plot differences in beta by iteration
#'
#' @param diffs
#'
#' @returns ggplot object
#' @export
#'
plot_differneces <- function(diffs){
  data.frame(diffs, id = 1:length(diffs)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = id, y = diffs)) +
    ggplot2::labs(title = "Differences",
                  x = "Iteration",
                  y = "Difference")
}



#' Plot the change in betas through the loop
#'
#' @param betas list of betas
#'
#' @returns ggplot object
#' @export
#'
beta_path <- function(betas, K, B, Z, time){
  if (missing(B)) {
    B <- get_B(time = time,
               order = 3,
               nknots = 3)
  }
  if (missing(Z)){
    Z <- matrix(1, nrow = length(time))
  }

  yhats <- purrr::map(betas, ~estimate_y(beta = .x,
                                         B = B,
                                         Z = Z,
                                         K = K)) %>%
    purrr::imap_dfr(~data.frame(time = time,
                                run = .y, .x))
  colnames(yhats) <- c("time", "run",
                       paste0("Taxa.", 1:K))

  yhats <- yhats %>%
    tidyr::pivot_longer(-c(time, run)) %>%
    dplyr::mutate(name = factor(name,
                                levels = paste0("Taxa.", 1:K)))

  yhats %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = time,
                                    y = value,
                                    group = run,
                                    alpha = run),
                       linewidth = .5) +
    ggplot2::facet_wrap(~name) +
    ggplot2::labs(title = "Y hats progression new",
                  x = "Time",
                  y = "ra")

}




#' Title
#'
#' @param res Cleverly model object
#' @param K Number of responses
#' @param tau MCP hyper parameter.
#' @param psi Hyperparameter for clustering penalty (larger drives pairwise differences to zero)
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param theta ADMM hyper parameter.
#' @param max_admm_iter Max number of iterations for the ADMM loop
#' @param max_outer_iter Max number of iterations for the outer loop (Algorithm 1)
#'
#' @returns ggplot object
#' @export
plot_clusters <- function(res, K, tau, psi, gammas, theta, max_admm_iter, max_outer_iter){
  cluster_df <- data.frame(
    K = factor(1:K, levels = 1:K),
    cluster = factor(res$clusters$membership))

  y_hat <- res$y_hat
  # plot clusters
  plot <- y_hat %>%
    dplyr::mutate(response = factor(response, levels = 1:K)) %>%
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
    ggplot2::labs(subtitle = paste0("psi:",psi,
                                    ", tau:",round(tau,2),
                                    ", theta:",theta,
                                    ", gamma:", gammas[1],",", gammas[2],
                                    ", admm:",max_admm_iter,
                                    ", outer:",max_outer_iter))

  return(plot)

}

#' Plot clusters based on a yhat data frame
#'
#' yhat should contain yhat, y, Z, and time
#'
#' @param yhat Matrix that contains yhat, y, Z, and time
#' @param chosen_cluster Cluster object
#' @param K Number of responses
#'
#' @returns ggplot object
#' @export
plot_clusters_yhat <- function(yhat, chosen_cluster, K = 12){

  cluster_df <- data.frame(
    K = factor(1:K, levels = 1:K),
    cluster = factor(chosen_cluster$membership))

  # plot clusters
  plot <- yhat %>%
    dplyr::mutate(response = factor(response, levels = 1:K)) %>%
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
      values = c(viridis::viridis(length(unique(cluster_df$cluster))), "grey50"),
      name = "Cluster"
    )

  return(plot)

}

#' plot_initial_fit
#'
#' @param res Cleverly model object
#' @param K Number of responses
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#'
#' @returns ggplot object
#' @export
plot_initial_fit <- function(res, K, gammas){
  # Initial fit:
  y_hat_init <- res$y_hat_init
  # plot clusters
  plot <- y_hat_init %>%
    dplyr::mutate(response = factor(response, levels = 1:K)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time)) +
    ggplot2::geom_point(ggplot2::aes(y = y,
                                     color = factor(Z),
                                     shape = factor(Z)),
                        size = .6, alpha = .8) +
    ggplot2::geom_line(ggplot2::aes(y = yhat,
                                    color = factor(Z),
                                    group = factor(Z)),
                       linewidth = 1) +
    ggplot2::facet_wrap(~response) +
    ggplot2::labs(subtitle   = paste0("gammas:",gammas))

  return(plot)
}


#' plot_cluster_path
#'
#' @param res Cleverly model object
#' @param psi Hyperparameter for clustering penalty (larger drives pairwise differences to zero)
#' @param tau MCP hyper parameter.
#' @param theta ADMM hyper parameter.
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param max_admm_iter Max number of iterations for the ADMM loop
#' @param max_outer_iter Max number of iterations for the outer loop (Algorithm 1)
#' @param duration
#'
#' @returns ggplot object
#' @export
plot_cluster_path <- function(res, psi, tau, theta, gammas, max_admm_iter, max_outer_iter, duration){
  # Cluster progress:
  cluster_track <- res$cluster_list
  plot <- purrr::imap_dfr(cluster_track,
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
                                   ", gamma:", gammas[1],",", gammas[2],
                                   ",admm_iter:",max_admm_iter,
                                   "outer_iter:",max_outer_iter,
                                   "time:", duration),
                  color = "Outer iteration")

  return(plot)
}


#' plot_alg2_convergence
#'
#' @param res Cleverly model object
#'
#' @returns ggplot object
#' @export
plot_alg2_convergence <- function(res){
  # alg2 convergence
  plot <- purrr::imap_dfr(res$alg_2_beta_diff,
                  ~data.frame(cluster = unlist(.x),
                              run = .y)) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(x = row, y = cluster, color = factor(run))) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Overall Iteration",
                  y = "ss beta",
                  title = "Algorithm 2 beta diffs",
                  color = "Outer iteration")

  return(plot)
}


#' plot_d_convergence
#'
#' @param res Cleverly model object
#'
#' @returns ggplot object
#' @export
plot_d_convergence <- function(res){
  plot <- purrr::imap_dfr(res$d_list,
                  ~data.frame(cluster = unlist(.x),
                              run = .y)) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(x = row, y = cluster, color = factor(run))) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Overall Iteration",
                  y = "d ",
                  title = "ADMM: d (differnce in v between s and s+1)",
                  color = "Outer iteration")
  return(plot)
}

#' plot_r_convergence
#'
#' @param res Cleverly model object
#'
#' @returns ggplot
#' @export
plot_r_convergence <- function(res){
  plot <- purrr::imap_dfr(res$r_list,
                  ~data.frame(cluster = unlist(.x),
                              run = .y)) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(x = row, y = cluster, color = factor(run))) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Overall Iteration",
                  y = "r",
                  title = "ADMM: r = beta_k - beta_k' - v_kappa",
                  color = "Outer iteration")

  return(plot)
}
