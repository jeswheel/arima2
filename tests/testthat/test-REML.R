test_that("arima: ML works", {
  set.seed(123)

  tol <- 1e-6

  x1 <- 10 + stats::arima.sim(model = list(ar = c(0.3)), n = 100)
  x2 <- -11 + stats::arima.sim(model = list(ar = c(0.3, -0.1), ma = c(0.8)), n = 100)

  expect_no_condition(arima(x1, order = c(1, 0, 0), method = "REML"))
  expect_no_condition(arima(x2, order = c(2, 0, 1), method = "REML"))
})

test_that("REML conditions work as expected", {
  x1 <- stats::arima.sim(model = list(ar = c(0.3)), n = 100)

  expect_warning(
    arima(x1, order = c(1, 0, 0), include.mean = FALSE, method = 'REML'),
    "Using 'method = \"CSS-ML\"' because REML is equivalent to ML when no intercept is present."
  )

})

test_that("Kalman REML objective matches exact dense matrix formula to machine precision", {
  y <- as.numeric(datasets::lh) # built-in deterministic series (n = 48)
  n <- length(y)
  phi <- 0.6
  theta <- -0.3

  # Exact dense matrix calculation in base R
  gamma0 <- (1 + 2 * phi * theta + theta^2) / (1 - phi^2)
  acf_vals <- stats::ARMAacf(ar = phi, ma = theta, lag.max = n - 1L)
  Sigma_unit <- stats::toeplitz(acf_vals) * gamma0 # Cov(y) / sigma_eps^2

  ones <- rep(1, n)
  S_11 <- as.numeric(crossprod(ones, solve(Sigma_unit, ones)))
  S_X1 <- as.numeric(crossprod(ones, solve(Sigma_unit, y)))
  mu_hat <- S_X1 / S_11
  res_c <- y - mu_hat

  S_Y <- as.numeric(crossprod(res_c, solve(Sigma_unit, res_c)))
  sigma2_reml <- S_Y / (n - 1)
  det_Sigma_unit <- as.numeric(determinant(Sigma_unit, logarithm = TRUE)$modulus)

  exact_reml_loglik <- -0.5 * (
    (n - 1) * log(2 * pi * sigma2_reml) + det_Sigma_unit + log(S_11) + (n - 1)
  )

  # Evaluate arima2 at fixed (phi, theta)
  fit_fixed <- arima(
    y, order = c(1, 0, 1), method = "REML", fixed = c(phi, theta)
  )

  expect_equal(fit_fixed$loglik, exact_reml_loglik, tolerance = 1e-12)
  expect_equal(unname(fit_fixed$coef["intercept"]), mu_hat, tolerance = 1e-12)
  expect_equal(fit_fixed$sigma2, sigma2_reml, tolerance = 1e-12)
})
