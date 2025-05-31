#' Cleverly: Main Function of the Package
#'
#' Runs the Cleverly algorithm over a range of \eqn{\psi} values and selects the optimal one based on BIC.
#'
#' @param Y A data frame or matrix of count response variables. Each column should represent a response variable, and each row a subject-time observation. Must be ordered by time. Total number of rows is \eqn{M}. Make sure Y counts are not transformed to RA or rarified.
#' @param Z A data frame or matrix containing external variables. Must have \eqn{M} rows and \eqn{L} columns.
#' @param subject_ids A vector of subject identifiers (length \eqn{M}) or a column name/index if \code{Y} is a data frame.
#' @param time A numeric vector of time points (length \eqn{M}) or a column name/index if \code{Y} is a data frame.
#' @param lp Either a numeric index or name of the column in \code{Z} used for clustering. Use 0 to indicate baseline clustering.
#' @param response_type Type of response variable. Should be either \code{"Counts"} or \code{"Continuous"}.
#' @param cor_str Correlation structure to use. One of \code{"IND"}, \code{"CON"}, \code{"AR1"}, \code{"CON-d"}, or \code{"AR1-d"}.
#' @param gammas A numeric vector of length \eqn{L + 1} specifying penalties applied to the D matrix.
#' @param psi_min Minimum \eqn{\psi} value (or the only value, if \code{npsi = 1}).
#' @param psi_max Maximum \eqn{\psi} value (ignored if \code{npsi = 1}).
#' @param npsi Number of \eqn{\psi} values to evaluate. If greater than 1, an equally spaced sequence from \code{psi_min} to \code{psi_max} is generated and the optimal value is selected via BIC.
#' @param parralel Logical; whether to use parallelization. Defaults to \code{FALSE}. If \code{TRUE}, uses \code{future::plan(future::multisession)}.
#' @param nworkers Number of workers to use if \code{parralel = TRUE}.
#' @param tau MCP hyperparameter. Default is \code{8/100}.
#' @param theta ADMM hyperparameter. Default is \code{300}.
#' @param C Constant controlling the Hessian update threshold. Default is \code{100}.
#' @param d Order of the finite difference matrix.
#' @param nknots Number of knots in the B-spline basis.
#' @param order Order of the B-spline basis.
#' @param epsilon_b Tolerance for convergence of Algorithm 1.
#' @param epsilon_r Tolerance for ADMM convergence (residual).
#' @param epsilon_d Tolerance for ADMM convergence (dual).
#' @param epsilon_2 Tolerance for convergence of Algorithm 2.
#' @param run_min Minimum number of runs for stability.
#' @param max_outer_iter Maximum number of iterations for the outer loop (Algorithm 1).
#' @param max_2_iter Maximum number of iterations for Algorithm 2 (per outer iteration).
#' @param max_admm_iter Maximum number of iterations for the ADMM clustering step (Algorithm 3).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{clusters}{Final cluster assignments.}
#'   \item{possible_clusters}{Cluster assignments for all tested \eqn{\psi} values.}
#'   \item{chosen_psi}{Optimal \eqn{\psi} value selected via BIC.}
#'   \item{y_hat_init}{Initial fitted values.}
#'   \item{y_hat}{Final fitted values.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- cleverly(
#'   Y = my_data[, c("response1", "response2")],
#'   Z = my_covariates,
#'   subject_ids = my_data$subject_id,
#'   time = my_data$time,
#'   lp = 1,
#'   response_type = "Counts",
#'   cor_str = "AR1",
#'   gammas = rep(1, ncol(my_covariates) + 1),
#'   psi_min = 0.1,
#'   psi_max = 1,
#'   npsi = 5,
#'   parralel = TRUE,
#'   nworkers = 4
#' )
#' }

cleverly <- function(Y,
                     Z,
                     subject_ids,
                     time,
                     lp = 0,
                     response_type = "counts",
                     cor_str = "IND",
                     gammas,
                     # optimize psi by BIC
                     psi_min = 500,
                     psi_max = 1500,
                     npsi = 10,
                     # parrallelilization options:
                     parralel = FALSE,
                     nworkers = 2,
                     # Other hyper parameters
                     tau = 0.005,
                     theta = 500,
                     C = 100,
                     d = 2,
                     nknots = 3,
                     order = 3,
                     # Convergence criteria
                     epsilon_b = 1e-3,
                     epsilon_r = 1e-3,
                     epsilon_d = 1e-3,
                     epsilon_2 = 1e-3,
                     # Run maximums
                     run_min = 3,
                     max_outer_iter = 30,
                     max_admm_iter = 100,
                     max_2_iter = 100,
                     BIC_type = "refit") {


# Hyperparameter checks ---------------------------------------------------



  if (tau * theta - 1 == 0) {
    stop("Tau and theta will cause MCP to be INF")
  }
  if (tau <= (1/theta)) {
    stop("No closed form solution for v with this tau theta combo")
  }



# Formating Y, id, time ---------------------------------------------------


  # Format Y
  Y_user <- Y

  id_quo <- rlang::enquo(subject_ids)
  time_quo <- rlang::enquo(time)

  if (is.data.frame(Y)) {
    # For ID:
    if (rlang::is_symbol(rlang::get_expr(id_quo)) &&
        rlang::as_name(id_quo) %in% colnames(Y)) {
      subject_ids <- dplyr::pull(Y, !!id_quo) # Extract the column
      Y <- dplyr::select(Y, -!!id_quo) # Remove the column
    } else if (length(rlang::eval_tidy(subject_ids)) == nrow(Y)) {
      subject_ids <- rlang::eval_tidy(id_quo) # Assume it's an external vector
    } else {
      stop("Invalid subject: must be a column name in Y or a vector of length nrow(Y)")
    }

    # For time:
    if (rlang::is_symbol(rlang::get_expr(time_quo)) &&
        rlang::quo_name(time_quo) %in% colnames(Y)) {
      time <- dplyr::pull(Y, !!time_quo) # Extract the column
      Y <- dplyr::select(Y, -!!time_quo) # Remove the column
    } else if (length(rlang::eval_tidy(time_quo)) == nrow(Y)) {
      time <- rlang::eval_tidy(time_quo) # Assume it's an external vector
    } else {
      stop("Invalid time: must be a column name in Y or a vector of length nrow(Y)")
    }
    # Convert Y into a matrix
    Y <- as.matrix(Y)
    dimnames(Y) <- NULL # To make the tests pass
  } else if (is.matrix(Y)) {
    if (length(subject_ids) != nrow(Y) || length(time) != nrow(Y)) {
      stop("For matrices, ID and Time must be external vectors of length nrow(Y)")
    }
  } else {
    stop("Y must be either a data frame or a matrix")
  }



# Format Z, mi ------------------------------------------------------------
  # Y_list <- get_Y_wrapper(Y = Y_user,
  #                         subject_ids = subject_ids,
  #                         time = time)
  # #Y <- Y_list$Y
  # # After this check time should be a numeric, and
  # time <- Y_list$time
  # subject_ids <- Y_list$subject_id_values
  # Get is and js
  # Time needs to be changed to indexes

  # Get m_is
  # Currently a data frame with subject id and mi
  mi_vec <- get_mi_vec(Y_user, subject_ids, time)$mi

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
  # Check gamma vector is of correct dimension
  if (length(gammas) != ncol(Z)) {
    stop("Length of gammas must be equal to the number of columns in Z + 1")
  }

# Run algorithm 1 ---------------------------------------------------------
  # Run algorithm 1 using correct formatting
  # Eventually this would change to a different algorithm for different response types
  if (response_type == "counts") {

    psis <- seq(psi_min, psi_max, length.out = npsi)

    if (parralel) {
      future::plan(future::multisession, workers = nworkers)
      res_list <- furrr::future_map(psis, ~algorithm1(Y = Y,
                                                      Z = Z,
                                                      time = time,
                                                      mi_vec = mi_vec,
                                                      lp = lp,
                                                      gammas = gammas,
                                                      psi = ..1,
                                                      tau = tau,
                                                      theta = theta,
                                                      C = C,
                                                      d = d,
                                                      run_min = run_min,
                                                      nknots = nknots,
                                                      order = order,
                                                      epsilon_b = epsilon_b,
                                                      epsilon_r = epsilon_r,
                                                      epsilon_d = epsilon_d,
                                                      max_outer_iter = max_outer_iter,
                                                      max_admm_iter = max_admm_iter,
                                                      max_2_iter = max_2_iter,
                                                      epsilon_2 = epsilon_2,
                                                      cor_str = cor_str))
      future::plan(future::sequential) # Reset to sequential after parallel processing
    } else {
      res_list <- list()
      for (p in 1:length(psis)) {
        psi <- psis[p]

        res_list[[p]] <- algorithm1(Y = Y,
                                    Z = Z,
                                    time = time,
                                    mi_vec = mi_vec,
                                    lp = lp,
                                    gammas = gammas,
                                    psi = psi,
                                    tau = tau,
                                    theta = theta,
                                    C = C,
                                    d = d,
                                    run_min = run_min,
                                    nknots = nknots,
                                    order = order,
                                    epsilon_b = epsilon_b,
                                    epsilon_r = epsilon_r,
                                    epsilon_d = epsilon_d,
                                    max_outer_iter = max_outer_iter,
                                    max_admm_iter = max_admm_iter,
                                    max_2_iter = max_2_iter,
                                    epsilon_2 = epsilon_2,
                                    cor_str = cor_str)
      }
    }

    # Calculate both unrefit (BIC) and refit (BIC_ra_group)
    BIC_list <- purrr::map(res_list, ~.x$BIC)
    BIC_ra_group <- purrr::map(res_list, ~.x$BIC_ra_group)

    # But just use the one specified
    if (BIC_type == "refit") {
      BICs <- purrr::map_dbl(res_list, ~.x$BIC_ra_group$BIC)
    } else {
      BICs <- purrr::map_dbl(res_list, ~.x$BIC$BIC)
    }
    # Which psi result had the minimum BIC
    best <- which.min(BICs)
    result <- res_list[[best]]
    result$all_clusters_psi <- purrr::map(res_list, ~.x$clusters)

    # What was the chosen cluster?
    clusters <- purrr::map(res_list, ~.x$clusters)
    print(paste0("psi:", psis,", cluster:", clusters))
    print(paste0("chosen psi: ", psis[best], ", cluster", clusters[[best]]))

  } else {
    stop("Invalid response type or type not yet implemented.")
  }

  # Return what is inputted and the result
  # Change eventually
  return(list(clusters = result$clusters,
              y_hat = result$y_hat,
              y_hat_init = result$y_hat_init,
              y_hat_baseline = result$y_hat_baseline,
              y_hat_lp_group = result$y_hat_lp_group,
              BIC = BIC_list,
              BIC_ra_group = BIC_ra_group,
              error = result$error,
              rho = result$rho,
              phi = result$phi,
              psi = psis[best],
              v = result$v,
              beta = result$beta,
              possible_clusters = clusters,
              s = result$s))



}





# Helpers: ----------------------------------------------------------------


#' Format Z
#'
#' Add a column of 1s to Z if it doesn't already exist
#'
#' @param Z A matrix or data frame with columns of external variables for each subject/time
#'
#' @returns A matrix with a column of 1s representing L = 0, and values for the other external variables
#' @export
format_Z <- function(Z) {
  if (is.data.frame(Z) | is.matrix(Z)) {
    M <- nrow(Z)
    if (!identical(Z[, 1], rep(1, M))) {
      Z <- cbind(1, Z)
    }
    Z <- as.matrix(Z)
  } else if (is.vector(Z)) {
    if (!all(Z == 1)) {
      Z <- cbind(1, Z)
    }
  }
  return(Z)
}



#' Title
#'
#' @param lp either a numeric index of which external variable to cluster on, or the name of the column of Z that contains the clustering variable. Specify numeric 0 to cluster via baseline.
#' @param Z Matrix that starts with a column of 1s. Of dimension M x (L + 1) that contains the external variable values for each subject/time and is 1 for l = 0. In the case that there are no external variables this is a matrix with one column of 1s.
#'
#' @returns lp index
#' @export
format_lp <- function(lp, Z) {
  return(lp)
}

