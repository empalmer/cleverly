

# Algorithm 3 -------------------------------------------------------------


# Get v, beta, lambda -----------------------------------------------------
test_that("update_v returns correct structure and dimensions", {
  set.seed(1)

  # Toy inputs
  K <- 3  # number of responses
  P <- 2  # number of spline coefficients
  L <- 2  # number of covariates (excluding intercept)
  lp <- 1 # use 0-based index per function doc

  beta <- matrix(rnorm(P * K * (L + 1)), nrow = P * K, ncol = L + 1)
  Kappa <- matrix(c(1, 2,
                    1, 3,
                    2, 3), ncol = 2, byrow = TRUE) # all pairwise combos

  lambda <- rnorm(nrow(Kappa) * P)
  gammas <- rep(0.1, L + 1)
  tau <- 0.5
  theta <- 1.0
  psi <- 0.3

  result <- update_v(beta, lp, lambda, Kappa, P, gammas, tau, theta, psi)

  # Basic structure checks
  expect_type(result, "list")
  expect_named(result, c("v", "u_list"))
  expect_true(is.numeric(result$v))
  expect_true(is.list(result$u_list))

  # Length checks
  expect_length(result$v, nrow(Kappa) * P)
  expect_length(result$u_list, nrow(Kappa))
  for (u in result$u_list) {
    expect_length(u, P)
  }
})

test_that("update_v handles zero norm_u case correctly", {
  # Construct beta so pairwise beta difference and lambda are zero
  P <- 2
  K <- 2
  L <- 1
  lp <- 0

  beta_val <- rep(0.5, P)
  beta <- matrix(0, nrow = P * K, ncol = L + 1)
  beta[1:P, lp + 1] <- beta_val
  beta[(P + 1):(2 * P), lp + 1] <- beta_val

  Kappa <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  lambda <- rep(0, P)
  gammas <- rep(0.1, L + 1)
  tau <- 0.5
  theta <- 1
  psi <- 0.3

  result <- update_v(beta, lp, lambda, Kappa, P, gammas, tau, theta, psi)

  # When norm_u == 0, v should equal u, which is 0
  expect_equal(result$v, rep(0, P))
})

