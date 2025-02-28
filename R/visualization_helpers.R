#' Visualize just the Ys.
#'
#' @param Y Matrix of counts. Each response should be a separate column (K). Each row should be a separate subject/time combination. There should be M total rows.
#' @param is numeric of length M of individual ids
#' @param time vector of time values for each subject/time
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#'
#' @returns ggplot object
#' @export
visualize_ys <- function(Y, is, time, beta){

  data.frame(time = time, id = factor(is), Y) %>%
    tidyr::pivot_longer(cols = -c(time, id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = value, color = id)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~name)
}


#' Plot (curve only) for a given beta
#'
#' @param beta matrix of beta (or beta hat) of dimension (P*K) x L
#' @param B B spline basis matrix of dimension (N x P)
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#' @param K Number of responses
#' @param time vector of time values for each subject/time
#'
#' @returns ggplot object
#' @export
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




