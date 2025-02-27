visualize_ys <- function(Y, is, time, beta){

  data.frame(time = time, id = factor(is), Y) %>%
    tidyr::pivot_longer(cols = -c(time, id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = value, color = id)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~name)
}


#' Plot (curve only) for a given beta
#'
#' @param beta
#' @param B
#' @param Z
#' @param K
#' @param time
#'
#' @returns
#' @export
#'
#' @examples
visualize_curve <- function(beta, B, Z, K, time){
  yhat <- estimate_y(beta = beta,
                      B = B,
                      Z = Z,
                      K = K)
  colnames(yhat) <- paste0("Taxa.", 1:K)

  yhat <- data.frame(time = time, yhat) %>%
    tidyr::pivot_longer(-c(time)) %>%
    dplyr::mutate(name = factor(name,
                                levels = paste0("Taxa.", 1:K)))

  plot <- yhat %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = time,
                                    y = value),
                       linewidth = .5) +
    ggplot2::facet_wrap(~name)

  return(plot)

}




beta_path <- function(betas){

  yhats <- purrr::map(betas, ~estimate_y(beta = .x,
                                         B = B,
                                         Z = Z,
                                         K = K)) %>%
    purrr::imap_dfr(~data.frame(time = time, run = .y, .x))
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
