test_that("GEEwrapper works", {
  result <- GEEclevR(Y,
                     Z,
                     subject_ids,
                     time,
                     lp,
                     gamma,
                     d = 3,
                     nknots = 3,
                     order = 3,
                     tol = 1e6,
                     smax = 100)
})
