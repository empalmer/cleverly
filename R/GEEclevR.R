#' GEEclevR
#'
#' @param Y either a data frame or matrix. Each response should be a separate column. Each row should be a separate subject/time combination. There should be M total rows.
#' @param subject_ids either a vector of length(Y) or a column reference if Y is a data frame
#' @param time either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric
#' @param lp
#' @param Z
#' @param gamma
#' @param d
#' @param nknots
#' @param order
#' @param tol
#' @param smax
#'
#' @returns
#' @export
#'
#' @examples
GEEclevR <- function(Y,
                     Z,
                     subject_ids,
                     time,
                     lp,
                     gamma,
                     d = 3,
                     nknots = 3,
                     order = 3,
                     tol = 1e6,
                     smax = 100) {

  # Checks ----------
  # This should actually be moved to a wrapper function! probably called geeclevR
  Y_user <- Y
  Y_list <- get_Y_wrapper(Y_user, subject_ids, time)
  mis <- get_mis(Y_user, subject_ids, time)

  Y <- Y_list$Y
  # After this check time should be a numeric, and
  time <- Y_list$time
  subject_ids <- Y_list$subject_ids

  Z_user <- Z
  Z <- format_Z(Z_user)


}
