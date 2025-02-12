test_that("extract yij", {
  Y <- structure(c(16L, 59L, 42L, 46L, 17L, 39L, 11L, 5L, 31L, 63L,
        27L, 52L, 11L, 39L, 27L, 28L, 8L, 4L, 52L, 22L, 55L, 67L, 43L,
        10L, 11L, 4L, 48L, 8L, 61L, 43L, 8L, 32L, 34L, 0L, 29L, 3L, 19L,
        46L, 6L, 21L, 28L, 0L, 72L, 0L, 13L, 48L, 1L, 20L, 2L, 32L, 3L,
        3L, 6L, 53L, 5L, 41L, 0L, 9L, 0L, 17L), dim = c(15L, 4L))
  expect_equal(get_Yij0(i = 1, j = 1, Y = Y, mi = 3), 100)
})
