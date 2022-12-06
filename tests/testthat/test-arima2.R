test_that("arima2: ML works", {
  set.seed(1)

  x1 <- arima.sim(model = list(ar = c(0.3)), n = 100)
  x2 <- arima.sim(model = list(ar = c(0.3, -0.1), ma = c(0.8)), n = 100)
  x3 <- arima.sim(model = list(ar = c(0.3), ma = c(-0.5, 0.2)), n = 100)

  arma2_1 <- arima2(x1, order = c(1, 0, 0), method = "ML")
  arma2_2 <- arima2(x2, order = c(2, 0, 1), method = "ML")
  arma2_3 <- arima2(x3, order = c(3, 0, 2), method = "ML")

  arma_1 <- stats::arima(x1, order = c(1, 0, 0), method = "ML")
  arma_2 <- stats::arima(x2, order = c(2, 0, 1), method = "ML")
  arma_3 <- stats::arima(x3, order = c(3, 0, 2), method = "ML")

  expect_gte(arma2_1$loglik, arma_1$loglik)
  expect_gte(arma2_2$loglik, arma_2$loglik)
  expect_gte(arma2_3$loglik, arma_3$loglik)
})

test_that("arima2: CSS-ML works", {
  set.seed(2)

  x1 <- arima.sim(model = list(ar = c(0.3)), n = 100)
  x2 <- arima.sim(model = list(ar = c(0.3, -0.1), ma = c(0.8)), n = 100)
  x3 <- arima.sim(model = list(ar = c(0.3), ma = c(-0.5, 0.2)), n = 100)

  arma2_1 <- arima2(x1, order = c(1, 0, 0))
  arma2_2 <- arima2(x2, order = c(2, 0, 1))
  arma2_3 <- arima2(x3, order = c(3, 0, 2))

  arma_1 <- stats::arima(x1, order = c(1, 0, 0))
  arma_2 <- stats::arima(x2, order = c(2, 0, 1))
  arma_3 <- stats::arima(x3, order = c(3, 0, 2))

  expect_gte(arma2_1$loglik, arma_1$loglik)
  expect_gte(arma2_2$loglik, arma_2$loglik)
  expect_gte(arma2_3$loglik, arma_3$loglik)
})
