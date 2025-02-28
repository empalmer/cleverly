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
