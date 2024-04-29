test_that("Sampling ARMA coefficients Works", {

  set.seed(123)
  Mod_bounds = c(0.05, 0.95)
  min_inv_root_dist = 0.1
  N <- 1000

  N_ar <- 3
  N_ma <- 3
  N_ar_seas <- 2
  N_ma_seas <- 1

  coefs <- sample_ARMA_coef(
    order = c(N_ar, N_ma),
    seasonal = list(order = c(N_ar_seas, N_ma_seas), period = 2),
    n = N, Mod_bounds = Mod_bounds, min_inv_root_dist = min_inv_root_dist
  )

  expect_equal(dim(coefs), c(N, N_ar + N_ma + N_ar_seas + N_ma_seas))

  ar_roots <- as.data.frame(t(apply(coefs, 1, ARMApolyroots)))
  ma_roots <- as.data.frame(t(apply(coefs, 1, ARMApolyroots, type = "MA")))

  expect_lte(max(Mod(as.matrix(1/ar_roots))), Mod_bounds[2])
  expect_gte(min(Mod(as.matrix(1/ar_roots))), Mod_bounds[1])
  expect_lte(max(Mod(as.matrix(1/ma_roots))), Mod_bounds[2])
  expect_gte(min(Mod(as.matrix(1/ma_roots))), Mod_bounds[1])

  check_dist <- function(x) {
    ma_inv_roots <- 1 / ARMApolyroots(x, type = 'MA')
    ar_inv_roots <- 1 / ARMApolyroots(x, type = 'AR')

    min(Mod(outer(ar_inv_roots, ma_inv_roots, FUN = '-')))
  }

  expect_gte(min(apply(coefs, 1, check_dist)), min_inv_root_dist)
})
