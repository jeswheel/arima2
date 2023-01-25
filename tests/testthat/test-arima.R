test_that("arima: ML works", {
  set.seed(1)

  x1 <- stats::arima.sim(model = list(ar = c(0.3)), n = 100)
  x2 <- stats::arima.sim(model = list(ar = c(0.3, -0.1), ma = c(0.8)), n = 100)
  x3 <- stats::arima.sim(model = list(ar = c(0.3), ma = c(-0.5, 0.2)), n = 100)

  arma2_1 <- arima(x1, order = c(1, 0, 0), method = "ML")
  arma2_2 <- arima(x2, order = c(2, 0, 1), method = "ML")
  arma2_3 <- arima(x3, order = c(3, 0, 2), method = "ML")

  arma_1 <- stats::arima(x1, order = c(1, 0, 0), method = "ML", SSinit = "Rossignol2011")
  arma_2 <- stats::arima(x2, order = c(2, 0, 1), method = "ML", SSinit = "Rossignol2011")
  arma_3 <- stats::arima(x3, order = c(3, 0, 2), method = "ML", SSinit = "Rossignol2011")

  expect_gte(arma2_1$loglik, arma_1$loglik)
  expect_gte(arma2_2$loglik, arma_2$loglik)
  expect_gte(arma2_3$loglik, arma_3$loglik)
})

test_that("arima: CSS-ML works", {
  set.seed(2)

  x1 <- stats::arima.sim(model = list(ar = c(0.3)), n = 100)
  x2 <- stats::arima.sim(model = list(ar = c(0.3, -0.1), ma = c(0.8)), n = 100)
  x3 <- stats::arima.sim(model = list(ar = c(0.3), ma = c(-0.5, 0.2)), n = 100)

  arma2_1 <- arima(x1, order = c(1, 0, 0))
  arma2_2 <- arima(x2, order = c(2, 0, 1))
  arma2_3 <- arima(x3, order = c(3, 0, 2))

  arma_1 <- stats::arima(x1, order = c(1, 0, 0), SSinit = "Rossignol2011")
  arma_2 <- stats::arima(x2, order = c(2, 0, 1), SSinit = "Rossignol2011")
  arma_3 <- stats::arima(x3, order = c(3, 0, 2), SSinit = "Rossignol2011")

  expect_gte(arma2_1$loglik, arma_1$loglik)
  expect_gte(arma2_2$loglik, arma_2$loglik)
  expect_gte(arma2_3$loglik, arma_3$loglik)
})

test_that("arima: CSS works", {
  set.seed(3)

  x <- arima.sim(model = list(ar = c(0.3, -0.2), ma = c(0.1)), n = 100)
  arma <- stats::arima(x, order = c(2, 0, 1), method = "CSS")
  arma2 <- arima(x, order = c(2, 0, 1), method = "CSS", SSinit = "Gardner1980")

  expect_equal(arma$coef, arma2$coef)
  expect_equal(arma$var.coef, arma2$var.coef)
  expect_equal(arma$loglik, arma2$loglik)
  expect_equal(arma$residuals, arma2$residuals)
})

test_that("arima: max_iters = 0 works", {
  set.seed(4)

  x <- arima.sim(model = list(ar = c(0.3, -0.2), ma = c(0.1)), n = 100)
  arma <- stats::arima(x, order = c(2, 0, 1))
  arma2 <- arima(x, order = c(2, 0, 1), max_iters = 1)

  expect_equal(arma$coef, arma2$coef)
  expect_equal(arma$var.coef, arma2$var.coef)
  expect_equal(arma$loglik, arma2$loglik)
  expect_equal(arma$residuals, arma2$residuals)
})
