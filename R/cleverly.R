#' cleverly
#'
#' @param Y either a data frame or matrix of numeric response variables. Each response should be a separate column. Each row should be a separate subject/time combination. There should be M total rows. Must be ordered in time.
#' @param Z Matrix or data frame containing a column for each external variable. There should be M rows and L columns.
#' @param subject_ids either a vector of length(Y) or a column reference if Y is a data frame
#' @param time either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric
#' @param lp either a numeric index of which external variable to cluster on, or the name of the column of Z that contains the clustering variable. Specify numeric 0 to cluster via baseline.
#' @param response_type
#' @param gamma Vector of length (L + 1) of penalization hyper parameters
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param tol Tolerance for convergence
#' @param smax Maximum number of iterations
#'
#' @returns
#' @export
#'
#' @examples
cleverly <- function(Y,
                     Z,
                     subject_ids,
                     time,
                     lp = 0,
                     response_type = "counts",
                     gamma,
                     phi,
                     d = 3,
                     nknots = 3,
                     order = 3,
                     tol = 1e6,
                     smax = 100) {

  # Format Y
  Y_user <- Y
  Y_list <- get_Y_wrapper(Y_user, subject_ids, time)
  Y <- Y_list$Y
  # After this check time should be a numeric, and
  time <- Y_list$time
  subject_ids <- Y_list$subject_id_values
  # Get is and js
  # Time needs to be changed to indexes

  is <- subject_ids



  # Get m_is
  mis <- get_mis(Y_user, subject_ids, time)
  js <- sequence(mis$mi)

  # Format Z
  if (!missing(Z)) {
    Z_user <- Z
    Z <- format_Z(Z_user)
  } else {
    Z <- matrix(1, nrow = nrow(Y))
  }

  # Format l_p
  if (lp != 0) {
    lp <- format_lp(lp = lp, Z = Z)
  }

  # Run algorithm 1 using correct formatting
  # Eventually this would change to a different algorithm for different response types
  if (response_type == "counts") {
    result <- algorithm1(Y = Y,
                         Z = Z,
                         is = is,
                         time = time,
                         mis = mis,
                         lp = lp,
                         gamma = gamma,
                         phi = phi,
                         d = d,
                         nknots = nknots,
                         order = order,
                         tol = tol,
                         smax = smax)
  } else {
    stop("Invalid response type or type not yet implemented.")
  }

  # Return what is inputted and the result
  # Change eventually
  return(list(result = result,
              Y = Y,
              mis = mis,
              time = time,
              subject_ids = subject_ids,
              Z = Z,
              lp = lp))
}
