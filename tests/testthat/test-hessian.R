test_that("Hessian il", {
  sim <- base_sim()

  hessian_il <- get_Hessian_il(i = 1,
                               l = 2,
                               Y = sim$Y,
                               mi_vec = sim$mi_vec,
                               beta = sim$beta,
                               Z = sim$Z,
                               B = sim$B,
                               phi = .5)

  expect_equal(dim(hessian_il),
               c(sim$K*sim$P, sim$K*sim$P))


})





test_that("Hessian l", {
  sim <- base_sim()

  hessian_l <- get_Hessian_l(l = 2,
                             Y = sim$Y,
                             mi_vec = sim$mi_vec,
                             beta = sim$beta,
                             Z = sim$Z,
                             B = sim$B,
                             phi = .5,
                             C = 20)

  expect_equal(dim(hessian_l),
               c(sim$K*sim$P, sim$K*sim$P))

})


