test_that("Sampling ARMA coefficients Works", {

  set.seed(1233)
  Mod_bounds = c(0.05, 0.95)
  N <- 100
  coefs <- sample_ARMA_coef(
    order = c(3, 3),
    seasonal = list(order = c(2, 1), period = 2),
    n = N, Mod_bounds = Mod_bounds
  )

  ar_roots <- as.data.frame(t(apply(coefs, 1, ARMApolyroots)))
  ma_roots <- as.data.frame(t(apply(coefs, 1, ARMApolyroots, type = "MA")))

  expect_lte(max(Mod(as.matrix(1/ar_roots))), Mod_bounds[2])
  expect_gte(min(Mod(as.matrix(1/ar_roots))), Mod_bounds[1])
  expect_lte(max(Mod(as.matrix(1/ma_roots))), Mod_bounds[2])
  expect_gte(min(Mod(as.matrix(1/ma_roots))), Mod_bounds[1])

  min_diff <- Inf
  for (i in 1:N) {
    diff <- min(Mod(outer(1 / unlist(ar_roots[i, ]), 1 / unlist(ma_roots[i, ]), FUN = '-')))
    if (diff < min_diff) min_diff <- diff
  }

  expect_gte(min_diff, 0.05)
})
