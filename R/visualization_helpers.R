visualize_ys <- function(Y, is, time, beta){

  data.frame(time = time, id = factor(is), Y) %>%
    tidyr::pivot_longer(cols = -c(time, id)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = value, color = id)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~name)
}
