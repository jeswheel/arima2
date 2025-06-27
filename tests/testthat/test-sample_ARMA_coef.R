test_that("Unif-Root Sampling Works", {

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
    n = N, Mod_bounds = Mod_bounds, min_inv_root_dist = min_inv_root_dist,
    method = 'UnifRoots'
  )

  expect_true(all(apply(coefs[, paste0('ma', 1:3)], 1, arima2:::.maCheck)))
  expect_true(all(apply(coefs[, paste0('ar', 1:3)], 1, arima2:::.arCheck)))
  expect_true(all(abs(coefs[, paste0('ma_seas', 1)]) < 1))
  expect_true(all(apply(coefs[, paste0('ar_seas', 1:2)], 1, arima2:::.arCheck)))

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

test_that("DL Sampling Works", {

  set.seed(123)
  min_inv_root_dist = 0.1
  N <- 1000

  N_ar <- 3
  N_ma <- 3
  N_ar_seas <- 2
  N_ma_seas <- 1

  coefs <- sample_ARMA_coef(
    order = c(N_ar, N_ma),
    seasonal = list(order = c(N_ar_seas, N_ma_seas), period = 2),
    n = N, min_inv_root_dist = min_inv_root_dist,
    method = 'DL'
  )

  expect_true(all(apply(coefs[, paste0('ma', 1:3)], 1, arima2:::.maCheck)))
  expect_true(all(apply(coefs[, paste0('ar', 1:3)], 1, arima2:::.arCheck)))
  expect_true(all(abs(coefs[, paste0('ma_seas', 1)]) < 1))
  expect_true(all(apply(coefs[, paste0('ar_seas', 1:2)], 1, arima2:::.arCheck)))

  expect_equal(dim(coefs), c(N, N_ar + N_ma + N_ar_seas + N_ma_seas))

  ar_roots <- as.data.frame(t(apply(coefs, 1, ARMApolyroots)))
  ma_roots <- as.data.frame(t(apply(coefs, 1, ARMApolyroots, type = "MA")))

  expect_lte(max(Mod(as.matrix(1/ar_roots))), 1)
  expect_lte(max(Mod(as.matrix(1/ma_roots))), 1)

  check_dist <- function(x) {
    ma_inv_roots <- 1 / ARMApolyroots(x, type = 'MA')
    ar_inv_roots <- 1 / ARMApolyroots(x, type = 'AR')

    min(Mod(outer(ar_inv_roots, ma_inv_roots, FUN = '-')))
  }

  expect_gte(min(apply(coefs, 1, check_dist)), min_inv_root_dist)
})

test_that("DL Mod Bounds Warning", {
  expect_warning(
    sample_ARMA_coef(
      order = c(2, 2),
      n = 5, min_inv_root_dist = 0.1, Mod_bounds = c(0.05, 0.95),
      method = "DL"
    )
  )
})
