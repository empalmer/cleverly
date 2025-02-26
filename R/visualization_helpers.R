visualize_ys <- function(Y, is, time, beta){

  data.frame(time = time, id = factor(is), Y) %>%
    tidyr::pivot_longer(cols = -c(time, id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = value, color = id)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~name)
}


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
