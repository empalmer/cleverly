#' Cleverly: main function of the package
#'
#' Runs cleverly one psi for a range of psi values
#'
#' @param Y either a data frame or matrix of numeric response variables. Each response should be a separate column. Each row should be a separate subject/time combination. There should be M total rows. Must be ordered in time.
#' @param Z Matrix or data frame containing a column for each external variable. There should be M rows and L columns.
#' @param subject_ids either a vector of length(Y) or a column reference if Y is a data frame
#' @param time either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric
#' @param lp either a numeric index of which external variable to cluster on, or the name of the column of Z that contains the clustering variable. Specify numeric 0 to cluster via baseline.
#' @param response_type Counts or continuous response
#' @param cor_str Type of correlation structure (IND, CON, AR1, CON-d, AR1-d)
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param psi_min
#' @param psi_max
#' @param npsi
#' @param parralel
#' @param nworkers
#' @param tau MCP hyper parameter. Default is 8/100.
#' @param theta ADMM hyper parameter. Default is 300.
#' @param C Constant for determining the hessian change. Default is 100.
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param epsilon_b Tolerance for alg 1 convergence
#' @param epsilon_r Tolerance for ADMM convergence
#' @param epsilon_d Tolerance for ADMM convergence
#' @param epsilon_2 Tolerance for convergence of algorithm 2
#' @param run_min minimum number of runs
#' @param max_outer_iter Number of iterations for the outer loop (Algorithm 1)
#' @param max_2_iter Maximum number of iterations for algorithm 2 to run each loop (Algorithm 2)
#' @param max_admm_iter Number of iterations for the clustering step (Algorithm 3)
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
                     theta = 300,
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
                     max_2_iter = 100) {

  psis <- seq(psi_min, psi_max, length.out = npsi)

  if (parralel) {
    future::plan(future::multisession, workers = nworkers)
    res_list <- furrr::future_map(psis, ~cleverly_onepsi(Y = Y,
                                                         Z = Z,
                                                         subject_ids = subject_ids,
                                                         lp = 0,
                                                         time = time,
                                                         # Hyperparameters
                                                         gammas = gammas,
                                                         tau = tau,
                                                         theta = theta,
                                                         psi = ..1,
                                                         C = 100,
                                                         run_min = run_min,
                                                         # Iterations max
                                                         max_admm_iter = max_admm_iter,
                                                         max_outer_iter = max_outer_iter,
                                                         max_2_iter = max_2_iter,
                                                         # Convergence criteria
                                                         epsilon_r = .001,
                                                         epsilon_d = .05,
                                                         epsilon_b = .01,
                                                         epsilon_2 = .001,
                                                         cor_str = cor_str))
    future::plan(future::sequential) # Reset to sequential after parallel processing
  } else {
    res_list <- list()
    for (p in 1:length(psis)) {
      psi <- psis[p]

      res_list[[p]] <- cleverly_onepsi(Y = Y,
                                       Z = Z,
                                       subject_ids = subject_ids,
                                       time = time,
                                       lp = lp,
                                       response_type = response_type,
                                       cor_str = cor_str,
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
                                       epsilon_2 = epsilon_2)
    }
  }

  best <- which.min(purrr::map_dbl(res_list, ~.x$BIC))
  res <- res_list[[best]]
  res$all_clusters_psi <- purrr::map(res_list, ~.x$clusters)

  clusters <- purrr::map(res_list, ~.x$clusters)
  print(paste0("all clusters: psi:", psis,", cluster:", clusters))
  print(paste0("chosen psi cluster", purrr::map_dbl(res_list, ~.x$clusters$no)[best]))
  print(paste0("chosen psi", psis[best], ", cluster", clusters[[best]]))

  result <- list("cluster" = res$clusters,
                 "possible_clusters" = res$all_clusters_psi,
                 "chosen_psi" = psis[best],
                 "y_hat_init" = res$y_hat_init,
                 "y_hat_final" = res$y_hat)

  return(result)

}






#' cleverly one psi
#'
#' Overall user interface for this algorithm for a fixed single value of psi
#' Formats data and calls `algorithm1()` to run the algorithm
#'
#'
#' @param Y either a data frame or matrix of numeric response variables. Each response should be a separate column. Each row should be a separate subject/time combination. There should be M total rows. Must be ordered in time.
#' @param Z Matrix or data frame containing a column for each external variable. There should be M rows and L columns.
#' @param subject_ids either a vector of length(Y) or a column reference if Y is a data frame
#' @param time either a vector of length(Y) or a column reference if Y is a data frame. Must be numeric
#' @param lp either a numeric index of which external variable to cluster on, or the name of the column of Z that contains the clustering variable. Specify numeric 0 to cluster via baseline.
#' @param response_type Counts or continuous response
#' @param d Order for the difference matrix
#' @param nknots Number of knots for the B-spline basis
#' @param order Order of the B-spline basis
#' @param gammas Vector of dimension L + 1 for penalizing the D matrix
#' @param psi Hyperparameter for clustering penalty (larger drives pairwise differences to zero)
#' @param tau MCP hyper parameter. Default is 8/100.
#' @param theta ADMM hyper parameter. Default is 300.
#' @param C Constant for determining the hessian change. Default is 10.
#' @param max_outer_iter Number of iterations for the outer loop (Algorithm 1)
#' @param max_admm_iter Number of iterations for the clustering step (Algorithm 3)
#' @param epsilon_b Tolerance for alg 1 convergence
#' @param epsilon_r Tolerance for ADMM convergence
#' @param epsilon_d Tolerance for ADMM convergence
#' @param max_2_iter Maximum number of iterations for algorithm 2 to run each loop
#' @param epsilon_2 Tolerance for convergence of algorithm 2
#' @param cor_str Type of correlation structure (IND, CON, AR1)
#' @param run_min
#'
#' @returns
#' @export
cleverly_onepsi <- function(Y,
                            Z,
                            subject_ids,
                            time,
                            lp = 0,
                            response_type = "counts",
                            cor_str = "IND",
                            gammas,
                            psi,
                            tau = 8/100,
                            theta = 300,
                            C = 10,
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
                            max_2_iter = 100) {


  if (tau * theta - 1 == 0) {
    stop("Tau and theta will cause MCP to be INF")
  }
  if (tau <= (1/theta)) {
    stop("No closed form solution for v with this tau theta combo")
  }
  # Format Y
  Y_user <- Y

  id_quo <- rlang::enquo(subject_ids)
  time_quo <- rlang::enquo(time)

  if (is.data.frame(Y)) {
    # For ID:
    if (rlang::quo_name(id_quo) %in% colnames(Y)) {
      subject_ids <- dplyr::pull(Y, !!id_quo) # Extract the column
      Y <- dplyr::select(Y, -!!id_quo) # Remove the column
    } else if (length(subject_ids) == nrow(Y)) {
      subject_ids <- subject_ids # Assume it's an external vector
    } else {
      stop("Invalid subject: must be a column name in Y or a vector of length nrow(Y)")
    }

    # For time:
    if (rlang::quo_name(time_quo) %in% colnames(Y)) {
      time <- dplyr::pull(Y, !!time_quo) # Extract the column
      Y <- dplyr::select(Y, -!!time_quo) # Remove the column
    } else if (length(time) == nrow(Y)) {
      time <- time # Assume it's an external vector
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



  # Y_list <- get_Y_wrapper(Y = Y_user,
  #                         subject_ids = subject_ids,
  #                         time = time)
  # #Y <- Y_list$Y
  # # After this check time should be a numeric, and
  # time <- Y_list$time
  # subject_ids <- Y_list$subject_id_values
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

  if (length(gammas) != ncol(Z)) {
    stop("Length of gammas must be equal to the number of columns in Z + 1")
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
  } else {
    stop("Invalid response type or type not yet implemented.")
  }

  # Return what is inputted and the result
  # Change eventually
  return(list(clusters = result$clusters,
              y_hat = result$y_hat,
              y_hat_init = result$y_hat_init,
              # beta = result$beta,
              # beta_init = result$beta_init,
              # v = result$v,
              # ts = result$ts,
              # rs = result$rs,
              # u_list = result$u_list,
              # admm_diffs = result$admm_diffs,
              # admm_beta_list = result$admm_beta_list,
              # phis_list = result$phis_list,
              # admm_beta_list = result$admm_beta_list,
              # loop_list_beta = result$loop_list_beta,
              # alg1_diff = result$alg1_diff,
              # alg_2_beta_diff = result$alg_2_beta_diff,
              # cluster_list = result$cluster_list,
              # r_list = result$r_list,
              # d_list = result$d_list,
              BIC = result$BIC,
              error = result$error
              # Input values
              # B = result$B,
              # Y = Y,
              # mi_vec = mi_vec,
              # time = time,
              # subject_ids = subject_ids,
              # Z = Z,
              # lp = lp,
              # rho_cor = result$rho_cor
              ))
}
