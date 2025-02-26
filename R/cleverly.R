#' cleverly
#'
#' Overall user interface for this algorithm
#' Formats data and calls `algorithm1()` to run the algorithm
#'
#'
#'
#' @param Y either a data frame or matrix of numeric response variables. Each response should be a separate column. Each row should be a separate subject/time combination. There should be M total rows. Must be ordered in time.
#' @param Z Matrix or data frame containing a column for each external variable. There should be M rows and L columns.
#' @param subject_ids either a vector of length(Y) or a column reference if Y is a data frame
#' @param time either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric
#' @param lp either a numeric index of which external variable to cluster on, or the name of the column of Z that contains the clustering variable. Specify numeric 0 to cluster via baseline.
#' @param response_type
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param psi MCP hyper parameter
#' @param phi Dirichlet Multinomial over dispersion parameter. Should be changed to be estimated.
#' @param tau MCP hyper parameter. Default is 8/100.
#' @param theta MCP hyper parameter. Default is 300.
#' @param C Constant for determining the hessian change. Default is 10.
#' @param max_outer_iter Number of iterations for the outer loop (Algorithm 1)
#' @param max_admm_iter Number of iterations for the clustering step (Algorithm 3)
#' @param tol Tolerance for convergence for the outer loop AND admm loop
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
                     gammas,
                     psi,
                     phi,
                     tau = 8/100,
                     theta = 300,
                     C = 10,
                     d = 2,
                     nknots = 3,
                     order = 3,
                     tol = 1e-6,
                     max_outer_iter = 10,
                     max_admm_iter = 100) {



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
  # Currently a data frame with subject id and mi
  mi_vec <- get_mi_vec(Y_user, subject_ids, time)$mi
  js <- sequence(mi_vec)


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
                         time = time,
                         mi_vec = mi_vec,
                         lp = lp,
                         gammas = gammas,
                         psi = psi,
                         phi = phi,
                         tau = tau,
                         theta = theta,
                         C = C,
                         d = d,
                         nknots = nknots,
                         order = order,
                         tol = tol,
                         max_outer_iter = max_outer_iter,
                         max_admm_iter = max_admm_iter)
  } else {
    stop("Invalid response type or type not yet implemented.")
  }


  # Return what is inputted and the result
  # Change eventually
  return(list(result = result,
              Y = Y,
              mi_vec = mi_vec,
              time = time,
              subject_ids = subject_ids,
              Z = Z,
              lp = lp))
}
