test_that("Y visualization", {
  sim <- base_sim()

  Y <- visualize_ys(Y = sim$Y,
                    is = sim$is,
                    time = sim$time)
})
