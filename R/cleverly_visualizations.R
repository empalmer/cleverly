# Plots for cleverly ------------------------------------------------

#' Plot clusters - old version
#'
#' For binary EV
#'
#' @param res cleverly model object
#' @param response_names vector of names for each response
#' @param scales "free_y" or "fixed" to control scaling of the y axis
#' @param order "response" or "cluster" to have original ordering or cluster ordering of responses on the facet plots
#' @param nrow number of rows for the facet plot
#' @param EV_color color for the EV line
#' @param curve_type = baseline, refit_baseline, or slope
#' @param Y_counts need if curve_type = "slope
#' @param Z_name Character string of column name for Z. Only needed when there are more than one Z columns in the data.
#'
#' @returns ggplot object
#' @export
plot_clusters <- function(res,
                          response_names = NULL,
                          order = "response",
                          curve_type = "baseline",
                          nrow = 3,
                          scales = "fixed",
                          EV_color = "grey50",
                          Y_counts = NULL,
                          Z_name = "Z"){

  if (is.null(response_names)) {
    response_names <- paste0("Response ", 1:length(res$clusters$membership))
  }

  if (curve_type == "baseline") {
    Y <- res$y_hat %>%
      dplyr::select(-yhat) %>%
      dplyr::rename(yhat = y_hat_baseline)

  }
  if (curve_type == "refit_baseline") {
    Y <- res$y_hat %>%
      dplyr::select(-yhat) %>%
      dplyr::rename(yhat = y_hat_group)
  }
  if (curve_type == "slope") {
    Y <- res$y_hat
  }

  #Z <- res$y_hat$Z


  Z_sym <- rlang::sym(Z_name)
  if (!Z_name %in% colnames(Y)) {
    stop("Invalid Z_name. Please provide a valid column name for Z. Check the colnames using colnames(res$Z)")
  }

  Z <- dplyr::pull(Y, !!Z_name)


  binary_Z <- length(unique(Z)) <= 2


  Y <- Y %>%
    dplyr::mutate(Z = !!Z_sym)
  if (binary_Z) {
    Y <- Y %>% dplyr::mutate(Z = factor(!!Z_sym))
  }




  response_names <- paste0("Cluster ", res$clusters$membership, " - ", response_names)



  # baseline ----------------------------------------------------------------
  if (curve_type == "baseline" | curve_type == "refit_baseline" ) {
    # Figure out how to order colors, based on response ordering or cluster ordering
    if (order == "response") {
      values <- viridis::viridis(length(unique(res$clusters$membership)))
    } else if (order == "cluster") {
      values <- sample(viridis::viridis(length(res$clusters$membership)))
    } else {
      stop("order must be either 'response' or 'cluster'")
    }


    Y$Z_color <- rep(min(as.numeric(as.character(Z))), length(Y$Z))
    if (binary_Z) {
      Y$Z_color <- factor(Y$Z_color)
    }
    plot <- Y %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names)) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                        color = Z,
                                        shape = cluster),
                           size = 1, alpha = .8) +
      ggplot2::labs(color = "Z",
                    shape = "Cluster") +
      ggplot2::geom_line(ggplot2::aes(y = yhat),
                         color = "black",
                         linewidth = 1.6) +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = Z_color),
                         #color = "#F8766D",
                         linewidth = 1.5) +
      #ggnewscale::new_scale_color() +
      # ggplot2::geom_line(ggplot2::aes(y = yhat,
      #                                 color = cluster),
      #                    linewidth = 1.5) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow)
    # ggplot2::scale_color_manual(values = values,
    #                             name = "Cluster")
  }


  # slope, binary ------------------------------------------------------------
  else if (binary_Z ) {
    # if (order == "response") {
    #   values <- c(viridis::viridis(length(unique(cluster_key$cluster))), EV_color)
    # } else if (order == "cluster") {
    #   values <-  c(sample(viridis::viridis(length(unique(cluster_key$cluster)))), EV_color)
    # } else {
    #   stop("order must be either 'response' or 'cluster'")
    #   # values = c(RColorBrewer::brewer.pal(n = length(unique(cluster_df$cluster)),
    #   #                       name = "Set1"),
    #   #            "grey30"),
    #   #values = c(rcartocolor::carto_pal(length(unique(cluster_df$cluster)),
    #   #                                  "Safe"),
    #   #           "grey50"),
    # }
    # plot clusters

    plot <- Y %>%
      dplyr::mutate(clusterZ = ifelse(.data$Z == 0, "Baseline Curve", .data$cluster),
                    response = factor(.data$response,
                                      labels = response_names)) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                        color = factor(Z),
                                        shape = cluster),
                           size = 1, alpha = .8) +
      ggplot2::guides(color = ggplot2::guide_legend("Z"),
                      shape = ggplot2::guide_legend("Cluster")) +
      # ggnewscale::new_scale_color() +
      ggplot2::geom_line(ggplot2::aes(y = yhat, color = factor(Z), group = Z),
                         linewidth = 1.5) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow) +
      ggplot2::scale_color_discrete(
        #values = values,
        name = "Cluster",
        labels = function(x) stringr::str_wrap(x, width = 10)
      )
    # slope, continuous ------------------------------------------------------------
  } else if (!binary_Z ) {


    response_val <- res$y_hat$response[1]
    Z_orig <- Z
    Z <- res$y_hat %>%
      dplyr::filter(response == response_val) %>%
      dplyr::pull(!!Z_name)
    time <- res$y_hat %>%
      dplyr::filter(response == response_val) %>%
      dplyr::pull(time)
    y_baseline <- res$y_hat_baseline$yhat

    # Round Z into 4 distinct values.
    bins <- cut(Z, breaks = 4)
    breaks <- levels(cut(Z, breaks = 4, include.lowest = TRUE)) # string intervals
    numeric_breaks <- as.numeric(gsub("\\(|\\]|\\[|\\)", "", unlist(strsplit(breaks, ","))))
    break_mat <- matrix(numeric_breaks, ncol = 2, byrow = TRUE)

    #midpoints <- rowMeans(as.data.frame(break_mat)) %>% round(2)
    mins <- apply(as.data.frame(break_mat), 1, min) %>% round(2)

    # Map factor levels to midpoints
    Z_mid <- mins[as.integer(bins)]

    Z_rounded <- data.frame(rep(1, length(Z_mid)), Z_mid)

    #need original counts of y
    y_hat <- estimate_y(beta = res$beta,
                        B = res$B,
                        Z = as.matrix(Z_rounded),
                        K = length(response_names),
                        Y = Y_counts,
                        time = time,
                        baseline = F)
  cluster_key <- data.frame(
    response_names = factor(1:length(response_names),
                            levels = 1:length(response_names)),
    cluster = factor(res$clusters$membership))
    y_hat <- y_hat %>%
      dplyr::mutate(response = factor(.data$response)) %>%
      dplyr::left_join(cluster_key, by = c("response" = "response_names")) %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names))

    if (order == "cluster") {
      cluster_key <- cluster_key %>%
        dplyr::arrange(.data$cluster)
    }


    plot <- y_hat %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names),
                    Z_orig = Z_orig,
                    y_baseline = y_baseline) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                        color = Z_orig,
                                        shape = cluster),
                           alpha = .6) +
      ggplot2::labs(color = "Z") +
      #ggplot2::guides(color = ggplot2::guide_legend("Z")) +
      #ggnewscale::new_scale_color() +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = Z_mid,
                                      group = factor(Z_mid)),
                         linewidth = 1) +
      ggplot2::guides(color = "none") +
      ggplot2::labs(shape = "Cluster") +
      # ggplot2::geom_line(ggplot2::aes(y = y_baseline), color = "grey",
      #                    linewidth = 1) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow)

  }

  return(plot)

}


#' Plot fits and data for a single cluster
#'
#'
#' @param res cleverly result
#' @param response_names vector of response naems
#' @param cluster_val numeric which cluster to plot
#' @param curve_type "baseline" "refit" or "slope"
#' @param nrow number of rows for the facet plot to include
#' @param scales "fixed" or "free_y"
#' @param Y_counts needed if curve_type = "slope" and continuous Z
#' @param Z_name Character string of column name for Z. Only needed when there are more than one Z columns in the data.
#'
#' @returns ggplot object for one cluster
#' @export
plot_one_cluster <- function(res,
                             cluster_val,
                             response_names = NULL,
                             curve_type = "baseline",
                             nrow = 3,
                             scales = "fixed",
                             Y_counts = NULL,
                             Z_name = "Z"){

  if (is.null(response_names)) {
    response_names <- paste0("Response ", 1:length(res$clusters$membership))
  }

  if (curve_type == "baseline") {
    Y <- res$y_hat %>%
      dplyr::select(-yhat) %>%
      dplyr::rename(yhat = y_hat_baseline)
  }
  if (curve_type == "refit_baseline") {
    Y <- res$y_hat %>%
      dplyr::select(-yhat) %>%
      dplyr::rename(yhat = y_hat_group)
  }
  if (curve_type == "slope") {
    Y <- res$y_hat
  }


  Z_sym <- rlang::sym(Z_name)
  if (!Z_name %in% colnames(Y)) {
    stop("Invalid Z_col. Please provide a valid column name for Z. Check the colnames using colnames(res$Z)")
  }

  Z <- dplyr::pull(Y, !!Z_name)

  binary_Z <- length(unique(Z)) <= 2


  Y <- Y %>%
    dplyr::mutate(Z = !!Z_sym)
  if (binary_Z) {
    Y <- Y %>% dplyr::mutate(Z = factor(!!Z_sym))
  }



  response_names <- paste0("Cluster ",
                           res$clusters$membership,
                           " - ",
                           response_names)
  # baseline ----------------------------------------------------------------


  if (curve_type == "baseline" | curve_type == "refit_baseline" ) {
    # Figure out how to order colors, based on response ordering or cluster ordering
    # values <- viridis::viridis(length(unique(cluster_key$cluster)))

    Y$Z_color <- rep(min(as.numeric(as.character(Z))), length(Y$Z))
    if (binary_Z) {
      Y$Z_color <- factor(Y$Z_color)
    }
    plot <- Y %>%
      dplyr::filter(cluster == cluster_val) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y, color = Z),
                           size = 1, alpha = .8) +
      ggplot2::guides(color = ggplot2::guide_legend("Z")) +
      ggplot2::geom_line(ggplot2::aes(y = yhat),
                         color = "black",
                         linewidth = 1.6) +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = Z_color),
                         #color = "#F8766D",
                         linewidth = 1.5) +
      # ggnewscale::new_scale_color() +
      # ggplot2::geom_line(ggplot2::aes(y = yhat,
      #                                 color = cluster),
      #                    linewidth = 1.5) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow) +
      ggplot2::labs(title = paste0("Baseline cluster ", cluster_val),
                    color = "Cluster")
  }


  # slope, binary ------------------------------------------------------------
  else if (binary_Z) {
    # plot clusters
    plot <- Y %>%
      dplyr::mutate(clusterZ = ifelse(.data$Z == 0, "Baseline Curve", .data$cluster),
                    response = factor(.data$response, labels = response_names)) %>%
      dplyr::filter(cluster == cluster_val) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                        color = factor(Z)),
                           size = 1, alpha = .8) +
      ggplot2::guides(color = ggplot2::guide_legend("Z"),
                      shape = ggplot2::guide_legend("Z")) +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = factor(Z), group = Z),
                         linewidth = 1.5) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow) +
      ggplot2::scale_color_discrete(
        name = "Cluster",
        labels = function(x) stringr::str_wrap(x, width = 10)
      ) +
      ggplot2::ggtitle(paste0("Slope cluster ", cluster_val))
    # slope, continuous ------------------------------------------------------------
  } else if (!binary_Z ) {


    response_val <- res$y_hat$response[1]
    Z_orig <- Z
    Z <- res$y_hat %>%
      dplyr::filter(response == response_val) %>%
      dplyr::pull(!!Z_name)
    time <- res$y_hat %>%
      dplyr::filter(response == response_val) %>%
      dplyr::pull(time)
    y_baseline <- res$y_hat_baseline$yhat

    # Round Z into 4 distinct values.
    bins <- cut(Z, breaks = 3)
    breaks <- levels(cut(Z, breaks = 3, include.lowest = TRUE)) # string intervals
    numeric_breaks <- as.numeric(gsub("\\(|\\]|\\[|\\)", "", unlist(strsplit(breaks, ","))))
    break_mat <- matrix(numeric_breaks, ncol = 2, byrow = TRUE)

    #midpoints <- rowMeans(as.data.frame(break_mat)) %>% round(2)
    mins <- apply(as.data.frame(break_mat), 1, min) %>% round(2)
    # Map factor levels to midpoints
    Z_mid <- mins[as.integer(bins)]

    Z_rounded <- data.frame(rep(1, length(Z_mid)), Z_mid)

    #need original counts of y
    y_hat <- estimate_y(beta = res$beta,
                        B = res$B,
                        Z = as.matrix(Z_rounded),
                        K = length(response_names),
                        Y = Y_counts,
                        time = time,
                        baseline = F)

    cluster_key <- data.frame(
      response_names = factor(1:length(response_names),
                              levels = 1:length(response_names)),
      cluster = factor(res$clusters$membership))
    y_hat <- y_hat %>%
      dplyr::mutate(response = factor(.data$response)) %>%
      dplyr::left_join(cluster_key, by = c("response" = "response_names")) %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names))



    plot <- y_hat %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names),
                    Z_orig = Z_orig,
                    y_baseline = y_baseline) %>%
      dplyr::filter(cluster == cluster_val) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                        color = Z_orig),
                           alpha = .6) +
      ggplot2::labs(color = "Z") +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = Z_mid,
                                      group = factor(Z_mid)),
                         linewidth = 1) +
      ggplot2::labs(shape = "Cluster") +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow) +
      ggplot2::ggtitle(paste0("Slope cluster ", cluster_val))

  }

  return(plot)

}


# Plots for cleverly, simulation/real data versions ------------------------------------------------

#' Plot clusters - old version
#'
#' For binary EV
#'
#' @param res cleverly model object
#' @param response_names vector of names for each response
#' @param scales "free_y" or "fixed" to control scaling of the y axis
#' @param order "response" or "cluster" to have original ordering or cluster ordering of responses on the facet plots
#' @param nrow number of rows for the facet plot
#' @param EV_color color for the EV line
#' @param curve_type = baseline, refit_baseline, or slope
#' @param Y_counts need if curve_type = "slope
#'
#' @returns ggplot object
#' @export
plot_clusters_old <- function(res,
                          response_names = NULL,
                          order = "response",
                          curve_type = "baseline",
                          nrow = 3,
                          scales = "fixed",
                          EV_color = "grey50",
                          Y_counts = NULL){

  if (is.null(response_names)) {
    response_names <- paste0("Response ", 1:length(res$clusters$membership))
  }

  if (curve_type == "baseline") {
    Y <- res$y_hat_baseline
  }
  if (curve_type == "refit_baseline") {
    Y <- res$y_hat_lp_group
  }
  if (curve_type == "slope") {
    Y <- res$y_hat
  }

  Z <- res$y_hat$Z

  binary_Z <- length(unique(Z)) <= 2

  if (binary_Z) {
    Y <- Y %>% dplyr::mutate(Z = factor(.data$Z))
  }


  # Data formatting, adding cluster info:
  cluster_key <- data.frame(
    response_names = factor(1:length(response_names),
                            levels = 1:length(response_names)),
    cluster = factor(res$clusters$membership))
  Y <- Y %>%
    dplyr::mutate(response = factor(.data$response)) %>%
    dplyr::left_join(cluster_key, by = c("response" = "response_names")) %>%
    dplyr::mutate(response = factor(.data$response, labels = response_names))

  response_names <- paste0("Cluster ", cluster_key$cluster, " - ", response_names)


# baseline ----------------------------------------------------------------


  if (curve_type == "baseline" | curve_type == "refit_baseline" ) {
    # Figure out how to order colors, based on response ordering or cluster ordering
    if (order == "response") {
      values <- viridis::viridis(length(unique(cluster_key$cluster)))
    } else if (order == "cluster") {
      values <- sample(viridis::viridis(length(unique(cluster_key$cluster))))
    } else {
      stop("order must be either 'response' or 'cluster'")
    }


    Y$Z_color <- rep(min(as.numeric(as.character(Z))), length(Y$Z))
    if (binary_Z) {
      Y$Z_color <- factor(Y$Z_color)
    }
    plot <- Y %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names)) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                        color = Z,
                                        shape = cluster),
                          size = 1, alpha = .8) +
      ggplot2::labs(color = "Z",
                    shape = "Cluster") +
      ggplot2::geom_line(ggplot2::aes(y = yhat),
                         color = "black",
                         linewidth = 1.6) +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = Z_color),
                         #color = "#F8766D",
                         linewidth = 1.5) +
      #ggnewscale::new_scale_color() +
      # ggplot2::geom_line(ggplot2::aes(y = yhat,
      #                                 color = cluster),
      #                    linewidth = 1.5) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow)
      # ggplot2::scale_color_manual(values = values,
      #                             name = "Cluster")
  }


# slope, binary ------------------------------------------------------------
  else if (binary_Z ) {
    # if (order == "response") {
    #   values <- c(viridis::viridis(length(unique(cluster_key$cluster))), EV_color)
    # } else if (order == "cluster") {
    #   values <-  c(sample(viridis::viridis(length(unique(cluster_key$cluster)))), EV_color)
    # } else {
    #   stop("order must be either 'response' or 'cluster'")
    #   # values = c(RColorBrewer::brewer.pal(n = length(unique(cluster_df$cluster)),
    #   #                       name = "Set1"),
    #   #            "grey30"),
    #   #values = c(rcartocolor::carto_pal(length(unique(cluster_df$cluster)),
    #   #                                  "Safe"),
    #   #           "grey50"),
    # }
    # plot clusters

    plot <- Y %>%
      dplyr::mutate(clusterZ = ifelse(.data$Z == 0, "Baseline Curve", .data$cluster),
                    response = factor(.data$response,
                                      labels = response_names)) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                       color = factor(Z),
                                       shape = cluster),
                          size = 1, alpha = .8) +
      ggplot2::guides(color = ggplot2::guide_legend("Z"),
                      shape = ggplot2::guide_legend("Cluster")) +
      # ggnewscale::new_scale_color() +
      ggplot2::geom_line(ggplot2::aes(y = yhat, color = factor(Z), group = Z),
                         linewidth = 1.5) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow) +
      ggplot2::scale_color_discrete(
        #values = values,
        name = "Cluster",
        labels = function(x) stringr::str_wrap(x, width = 10)
      )
# slope, continuous ------------------------------------------------------------
  } else if (!binary_Z ) {

    response_val <- res$y_hat$response[1]
    Z_orig <- Z
    Z <- res$y_hat %>%
      dplyr::filter(response == response_val) %>%
      dplyr::pull(Z)
    time <- res$y_hat %>%
      dplyr::filter(response == response_val) %>%
      dplyr::pull(time)
    y_baseline <- res$y_hat_baseline$yhat

    # Round Z into 4 distinct values.
    bins <- cut(Z, breaks = 4)
    breaks <- levels(cut(Z, breaks = 4, include.lowest = TRUE)) # string intervals
    numeric_breaks <- as.numeric(gsub("\\(|\\]|\\[|\\)", "", unlist(strsplit(breaks, ","))))
    break_mat <- matrix(numeric_breaks, ncol = 2, byrow = TRUE)

    #midpoints <- rowMeans(as.data.frame(break_mat)) %>% round(2)
    mins <- apply(as.data.frame(break_mat), 1, min) %>% round(2)

    # Map factor levels to midpoints
    Z_mid <- mins[as.integer(bins)]

    Z_rounded <- data.frame(rep(1, length(Z_mid)), Z_mid)

    #need original counts of y
    y_hat <- estimate_y(beta = res$beta,
                        B = res$B,
                        Z = as.matrix(Z_rounded),
                        K = length(response_names),
                        Y = Y_counts,
                        time = time,
                        baseline = F)

    y_hat <- y_hat %>%
      dplyr::mutate(response = factor(.data$response)) %>%
      dplyr::left_join(cluster_key, by = c("response" = "response_names")) %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names))

    if (order == "cluster") {
      cluster_key <- cluster_key %>%
        dplyr::arrange(.data$cluster)
    }

    plot <- y_hat %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names),
                    Z_orig = Z_orig,
                    y_baseline = y_baseline) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                        color = Z_orig,
                                        shape = cluster),
                          alpha = .6) +
      ggplot2::labs(color = "Z") +
      #ggplot2::guides(color = ggplot2::guide_legend("Z")) +
      #ggnewscale::new_scale_color() +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = Z_mid,
                                      group = factor(Z_mid)),
                         linewidth = 1) +
      ggplot2::guides(color = "none") +
      ggplot2::labs(shape = "Cluster") +
      # ggplot2::geom_line(ggplot2::aes(y = y_baseline), color = "grey",
      #                    linewidth = 1) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow)

  }

  return(plot)

}

#' Plot fits and data for a single cluster
#'
#' @param res cleverly result
#' @param response_names vector of response naes
#' @param cluster_val numeric which cluster to plot
#' @param curve_type "baseline" "refit" or "slope"
#' @param nrow number of rows for the facet plot to include
#' @param scales "fixed" or "free_y"
#' @param Y_counts needed if curve_type = "slope" and continuous Z
#'
#' @returns ggplot object for one cluster
#' @export
plot_one_cluster_old <- function(res,
                             cluster_val,
                             response_names = NULL,
                             curve_type = "baseline",
                             nrow = 3,
                             scales = "fixed",
                             Y_counts = NULL){

  if (is.null(response_names)) {
    response_names <- paste0("Response ", 1:length(res$clusters$membership))
  }

  if (curve_type == "baseline") {
    Y <- res$y_hat_baseline
  }
  if (curve_type == "refit_baseline") {
    Y <- res$y_hat_lp_group
  }
  if (curve_type == "slope") {
    Y <- res$y_hat
  }

  Z <- Y$Z
  binary_Z <- length(unique(Z)) <= 2

  if (binary_Z ) {
    Y <- Y %>% dplyr::mutate(Z = factor(.data$Z))
  }


  # Data formatting, adding cluster info:
  cluster_key <- data.frame(
    response_names = factor(1:length(response_names),
                            levels = 1:length(response_names)),
    cluster = factor(res$clusters$membership))
  Y <- Y %>%
    dplyr::mutate(response = factor(.data$response)) %>%
    dplyr::left_join(cluster_key, by = c("response" = "response_names")) %>%
    dplyr::mutate(response = factor(.data$response, labels = response_names))

  response_names <- paste0("Cluster ",
                           cluster_key$cluster,
                           " - ",
                           response_names)
  # baseline ----------------------------------------------------------------


  if (curve_type == "baseline" | curve_type == "refit_baseline" ) {
    # Figure out how to order colors, based on response ordering or cluster ordering
    # values <- viridis::viridis(length(unique(cluster_key$cluster)))

    Y$Z_color <- rep(min(as.numeric(as.character(Z))), length(Y$Z))
    if (binary_Z) {
      Y$Z_color <- factor(Y$Z_color)
    }
    plot <- Y %>%
      dplyr::filter(cluster == cluster_val) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y, color = Z),
                           size = 1, alpha = .8) +
      ggplot2::guides(color = ggplot2::guide_legend("Z")) +
      ggplot2::geom_line(ggplot2::aes(y = yhat),
                         color = "black",
                         linewidth = 1.6) +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = Z_color),
                         #color = "#F8766D",
                         linewidth = 1.5) +
      # ggnewscale::new_scale_color() +
      # ggplot2::geom_line(ggplot2::aes(y = yhat,
      #                                 color = cluster),
      #                    linewidth = 1.5) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow) +
      ggplot2::labs(title = paste0("Baseline cluster ", cluster_val),
                    color = "Cluster")
  }


  # slope, binary ------------------------------------------------------------
  else if (binary_Z) {

    #values <- c(viridis::viridis(1), EV_color)

    # plot clusters
    plot <- Y %>%
      dplyr::mutate(clusterZ = ifelse(.data$Z == 0, "Baseline Curve", .data$cluster),
                    response = factor(.data$response, labels = response_names)) %>%
      dplyr::filter(cluster == cluster_val) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                       color = factor(Z)),
                          size = 1, alpha = .8) +
      ggplot2::guides(color = ggplot2::guide_legend("Z"),
                      shape = ggplot2::guide_legend("Z")) +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = factor(Z), group = Z),
                         linewidth = 1.5) +
      # ggnewscale::new_scale_color() +
      # ggplot2::geom_line(ggplot2::aes(y = yhat,
      #                                 color = clusterZ,
      #                                 group = factor(Z)),
      #                    linewidth = 1) +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow) +
      ggplot2::scale_color_discrete(
        #values = values,
        name = "Cluster",
        labels = function(x) stringr::str_wrap(x, width = 10)
      ) +
      ggplot2::ggtitle(paste0("Slope cluster ", cluster_val))
    # slope, continuous ------------------------------------------------------------
  } else if (!binary_Z ) {

    response_val <- res$y_hat$response[1]
    Z_orig <- res$y_hat$Z
    Z <- res$y_hat %>%
      dplyr::filter(response == response_val) %>%
      dplyr::pull(Z)
    time <- res$y_hat %>%
      dplyr::filter(response == response_val) %>%
      dplyr::pull(time)
    y_baseline <- res$y_hat_baseline$yhat

    # Round Z into 4 distinct values.
    bins <- cut(Z, breaks = 3)
    breaks <- levels(cut(Z, breaks = 3, include.lowest = TRUE)) # string intervals
    numeric_breaks <- as.numeric(gsub("\\(|\\]|\\[|\\)", "", unlist(strsplit(breaks, ","))))
    break_mat <- matrix(numeric_breaks, ncol = 2, byrow = TRUE)

    #midpoints <- rowMeans(as.data.frame(break_mat)) %>% round(2)
    mins <- apply(as.data.frame(break_mat), 1, min) %>% round(2)
    # Map factor levels to midpoints
    Z_mid <- mins[as.integer(bins)]

    Z_rounded <- data.frame(rep(1, length(Z_mid)), Z_mid)

    #need original counts of y
    y_hat <- estimate_y(beta = res$beta,
                        B = res$B,
                        Z = as.matrix(Z_rounded),
                        K = length(response_names),
                        Y = Y_counts,
                        time = time,
                        baseline = F)

    y_hat <- y_hat %>%
      dplyr::mutate(response = factor(.data$response)) %>%
      dplyr::left_join(cluster_key, by = c("response" = "response_names")) %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names))




    plot <- y_hat %>%
      dplyr::mutate(response = factor(.data$response, labels = response_names),
                    Z_orig = Z_orig,
                    y_baseline = y_baseline) %>%
      dplyr::filter(cluster == cluster_val) %>%
      ggplot2::ggplot(ggplot2::aes(x = time)) +
      ggplot2::geom_jitter(ggplot2::aes(y = y,
                                        color = Z_orig),
                           alpha = .6) +
      ggplot2::labs(color = "Z") +
      ggplot2::geom_line(ggplot2::aes(y = yhat,
                                      color = Z_mid,
                                      group = factor(Z_mid)),
                         linewidth = 1) +
      ggplot2::labs(shape = "Cluster") +
      ggplot2::facet_wrap(~response,
                          scales = scales,
                          nrow = nrow) +
      ggplot2::ggtitle(paste0("Slope cluster ", cluster_val))

  }

  return(plot)

}





#' Plot BIC values for a result.
#'
#' @param res cleverly object
#' @param BIC_type which version of BIC to plot
#' @param psis sequence of psis used
#'
#' @returns ggplot object
#' @export
plot_BIC <- function(res, BIC_type = "BIC", psis){

  if (BIC_type == "BIC") {
    data <- res$BIC
  }
  if (BIC_type == "BIC_ra_group") {
    data <- res$BIC_ra_group
  }
  if (BIC_type == "BIC_group") {
    data <- res$BIC_group
  }

  possible_clusters <- res$possible_clusters %>% purrr::map_dbl("no")
  bics <- purrr::map_dbl(data, "BIC")
  first_term <- purrr::map_dbl(data, "first_term")
  second_term <- purrr::map_dbl(data, "second_term")

  chosen <- which(psis == res$psi)

  df <- data.frame(possible_clusters,
                   bics,
                   first_term = first_term,
                   second_term)

  df$row_id <- seq_len(nrow(df))

  row_id <- seq_len(nrow(df))
  possible_clusters <- df$possible_clusters

  plot <- df %>%
    tidyr::pivot_longer(-c(possible_clusters, row_id), names_to = "term", values_to = "value") %>%
    ggplot2::ggplot(ggplot2::aes(x = row_id)) +
    ggplot2::geom_line(ggplot2::aes(y = value, color = term)) +
    ggplot2::scale_x_continuous(
      breaks = row_id,
      labels = paste0("C: ", possible_clusters, ", psi: ", round(psis, 2))
    ) +
    ggplot2::facet_wrap(~term, nrow = 3, scales = "free") +
    ggplot2::geom_vline(xintercept = chosen, linetype = "dashed", color = "red") +
    ggplot2::labs(x = "Possible Clusters (per row)", y = "BIC",
         title = paste0("BIC using: ", BIC_type))

  return(plot)
}





#' plot_initial_fit
#'
#' @param res Cleverly model object
#' @param response_names Vector of response names.
#' @param Z_name Character string of column name for Z. Default is "Z".
#'
#' @returns ggplot object
#' @export
plot_initial_fit <- function(res,
                             response_names = NULL,
                             Z_name = "Z"){

  if (is.null(response_names)) {
    response_names <- paste0("Response ", 1:length(res$clusters$membership))
  }

  K <- length(res$clusters$membership)
  # Initial fit:
  y_hat_init <- res$y_hat
  Z_sym <- rlang::sym(Z_name)


  if (!Z_name %in% colnames(y_hat_init)) {
    stop("Invalid Z_name. Please provide a valid column name for Z.")
  }

  binary_Z <- length(as.numeric(unique(y_hat_init[[Z_name]]))) <= 2

  if (binary_Z) {
    y_hat_init <- y_hat_init %>%
      dplyr::mutate(Z = factor(!!Z_sym))
  }

  # plot clusters
  plot <- y_hat_init %>%
    dplyr::mutate(response = factor(.data$response,
                                    levels = 1:K,
                                    labels = response_names)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time)) +
    ggplot2::geom_jitter(ggplot2::aes(y = .data$y,
                                     color = !!Z_sym),
                        size = .6, alpha = .8) +
    ggplot2::geom_line(ggplot2::aes(y = .data$y_hat_init,
                                    color = !!Z_sym,
                                    group = !!Z_sym),
                       linewidth = 1) +
    ggplot2::facet_wrap(~response) +
    ggplot2::labs(x = "Time", y = "Abundance", color = "Z")

  return(plot)
}


#' Plot an alluvial plot to compare two cluster results
#'
#' @param res1 cluster result 1
#' @param res2 cluster result 2
#' @param response_names vector of response names
#' @param res_names vector of length two to name the cluster results
#'
#' @returns ggplot object
#' @export
plot_cluster_differences <- function(res1,
                                     res2,
                                     response_names = NULL,
                                     res_names = c("A", "B")) {

  if (is.null(response_names)) {
    response_names <- paste0("Response ", 1:length(res1$clusters$membership))
  }
  # compare baseline clusters for high nutrition and body condition
  df <- data.frame(
    id = seq_along(res1$clusters$membership),
    response_names = response_names,
    res1 = as.factor(res1$clusters$membership),
    res2 = as.factor(res2$clusters$membership)
  )
  #
  #   ggplot2::ggplot(df, ggplot2::aes(axis1 = res1, axis2 = res2, y = 1)) +
  #     ggalluvial::geom_alluvium(ggplot2::aes(fill = res1),
  #                               width = 1/12, alpha = 0.6) +
  #     ggalluvial::geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  #     ggplot2::geom_text(ggplot2::aes(label = response_names),
  #                        stat = ggalluvial::StatAlluvium, size = 2) +
  #     ggplot2::scale_x_discrete(limits = res_names, expand = c(.05, .05)) +
  #     ggplot2::theme_minimal() +
  #     ggplot2::theme(axis.text.y = ggplot2::element_blank(),
  #                    axis.ticks.y = ggplot2::element_blank(),
  #                    axis.title.y = ggplot2::element_blank(),
  #                    legend.position = "none")

  df <- df[order(df$res1, df$res2, df$response_names), ]
  df$flow_id <- factor(df$response_names, levels = unique(df$response_names))

  df$response_names <- factor(df$response_names,
                              levels = sort(unique(df$response_names)))


  plot <- ggplot2::ggplot(df, ggplot2::aes(axis1 = res1,
                                           axis2 = res2,
                                           y = 1,
                                           group = response_names)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = res1),
                              width = 1/2,
                              alpha = 0.9) +
    ggalluvial::geom_stratum(width = 1/4,
                             fill = "grey90",
                             color = "black",
                             size = 1) +
    ggplot2::geom_text(stat = "alluvium",
                       ggplot2::aes(label = as.character(response_names)),
                       size = 2.5,
                       hjust = 0,
                       nudge_x = -0.12) +
    ggplot2::theme_void() +
    ggplot2::scale_x_discrete(
      limits = res_names,
      expand = c(.05, .05),
      position = "top"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(b = 2)),
      axis.title.x = ggplot2::element_text(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )

  return(plot)

}



# Plots for understanding algorithm  -------------------------------------------------------------------

#' Plot the final fit
#'
#' Plot the final fit of the model (relative abundance scale) with the data points (relative abundance scale)
#'
#' @param Y Input matrix Matrix of counts Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param res cleverly model output object
#'
#' @returns ggplot object
visualize_final_fit <- function(Y, res){
  K <- ncol(Y)

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
    dplyr::mutate(taxa = factor(.data$response,
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
#' @param diffs differences
#'
#' @returns ggplot object

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
#' @param K number of responses
#' @param B bspline
#' @param Z EV
#' @param time time
#'
#' @returns ggplot object
beta_path <- function(betas, K, B, Z, time){
  if (missing(B)) {
    B <- get_B(time = time,
               order = 3,
               nknots = 3)
  }
  if (missing(Z)) {
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
    tidyr::pivot_longer(-c(time, .data$run)) %>%
    dplyr::mutate(name = factor(.data$name,
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





#' Plot clusters based on a yhat data frame
#'
#' yhat should contain yhat, y, Z, and time
#'
#' @param yhat Matrix that contains yhat, y, Z, and time
#' @param chosen_cluster Cluster object
#' @param K Number of responses
#'
#' @returns ggplot object
plot_clusters_yhat <- function(yhat, chosen_cluster, K = 12){

  cluster_df <- data.frame(
    K = factor(1:K, levels = 1:K),
    cluster = factor(chosen_cluster$membership))

  # plot clusters
  plot <- yhat %>%
    dplyr::mutate(response = factor(.data$response, levels = 1:K)) %>%
    dplyr::left_join(cluster_df, by = c("response" = "K")) %>%
    dplyr::mutate(clusterZ = ifelse(.data$Z == 1, "EV", .data$cluster)) %>%
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




#' plot_cluster_path
#'
#' @param res Cleverly model object
#' @param psi Hyperparameter for clustering penalty (larger drives pairwise differences to zero)
#' @param tau MCP hyper parameter.
#' @param theta ADMM hyper parameter.
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param max_admm_iter Max number of iterations for the ADMM loop
#' @param max_outer_iter Max number of iterations for the outer loop (Algorithm 1)
#' @param duration duration
#'
#' @returns ggplot object
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
