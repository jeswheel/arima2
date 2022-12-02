#' Auto Arima 2
#' @inherit forecast::auto.arima
#' @param nrestart number of random restarts for arima2 function.
#' @export
auto.arima2 = function (y, nrestart = 10, d = NA, D = NA, max.p = 5, max.q = 5, max.P = 2,
          max.Q = 2, max.order = 5, max.d = 2, max.D = 1, start.p = 2,
          start.q = 2, start.P = 1, start.Q = 1, stationary = FALSE,
          seasonal = TRUE, ic = c("aicc", "aic", "bic"), stepwise = TRUE,
          nmodels = 94, trace = FALSE, approximation = (length(x) >
                                                          150 | stats::frequency(x) > 12), method = NULL, truncate = NULL,
          xreg = NULL, test = c("kpss", "adf", "pp"), test.args = list(),
          seasonal.test = c("seas", "ocsb", "hegy", "ch"), seasonal.test.args = list(),
          allowdrift = TRUE, allowmean = TRUE, lambda = NULL, biasadj = FALSE,
          parallel = FALSE, num.cores = 2, x = y, ...)
{

  ndiffs <- utils::getFromNamespace("ndiffs", "forecast")
  fitted.Arima <- utils::getFromNamespace("fitted.Arima", "forecast")
  UndoWhichModels <- utils::getFromNamespace("UndoWhichModels", "forecast")
  WhichModels <- utils::getFromNamespace("WhichModels", "forecast")
  if (stepwise && parallel) {
    warning("Parallel computer is only implemented when stepwise=FALSE, the model will be fit in serial.")
    parallel <- FALSE
  }
  if (trace && parallel) {
    message("Tracing model searching in parallel is not supported.")
    trace <- FALSE
  }
  series <- deparse(substitute(y))
  x <- stats::as.ts(x)
  if (NCOL(x) > 1) {
    stop("auto.arima2 can only handle univariate time series")
  }
  orig.x <- x
  missing <- is.na(x)
  firstnonmiss <- utils::head(which(!missing), 1)
  lastnonmiss <- utils::tail(which(!missing), 1)
  serieslength <- sum(!missing[firstnonmiss:lastnonmiss])
  x <- subset(x, start = firstnonmiss)
  if (!is.null(xreg)) {
    if (!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- as.matrix(xreg)
    xreg <- xreg[firstnonmiss:NROW(xreg), , drop = FALSE]
  }
  if (is.constant(x)) {
    if (all(is.na(x)))
      stop("All data are missing")
    if (allowmean) {
      fit <- arima2(x, order = c(0, 0, 0), fixed = mean(x,
                                                       na.rm = TRUE), nrestart=nrestart, ...)
    }
    else {
      fit <- arima2(x, order = c(0, 0, 0), include.mean = FALSE, nrestart=nrestart,
                   ...)
    }
    fit$x <- orig.x
    fit$series <- series
    fit$call <- match.call()
    fit$call$x <- data.frame(x = x)
    fit$constant <- TRUE
    return(fit)
  }
  ic <- match.arg(ic)
  test <- match.arg(test)
  seasonal.test <- match.arg(seasonal.test)
  if (seasonal) {
    m <- stats::frequency(x)
  }
  else {
    m <- 1
  }
  if (m < 1) {
    m <- 1
  }
  else {
    m <- round(m)
  }
  max.p <- min(max.p, floor(serieslength/3))
  max.q <- min(max.q, floor(serieslength/3))
  max.P <- min(max.P, floor(serieslength/3/m))
  max.Q <- min(max.Q, floor(serieslength/3/m))
  if (serieslength <= 3L) {
    ic <- "aic"
  }
  if (!is.null(lambda)) {
    x <- forecast::BoxCox(x, lambda)
    lambda <- attr(x, "lambda")
    attr(lambda, "biasadj") <- biasadj
  }
  if (!is.null(xreg)) {
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- if (ncol(xreg) == 1)
        "xreg"
      else paste("xreg", 1:ncol(xreg), sep = "")
    }
    xregg <- xreg
    xx <- x
    constant_columns <- apply(xregg, 2, is.constant)
    if (all(constant_columns)) {
      xregg <- NULL
    }
    else {
      if (any(constant_columns)) {
        xregg <- xregg[, -which(constant_columns), drop = FALSE]
      }
      sv <- svd(stats::na.omit(cbind(rep(1, NROW(xregg)), xregg)))$d
      if (min(sv)/sum(sv) < .Machine$double.eps) {
        stop("xreg is rank deficient")
      }
      j <- !is.na(x) & !is.na(rowSums(xregg))
      xx[j] <- stats::residuals(stats::lm(x ~ xregg))
    }
  }
  else {
    xx <- x
    xregg <- NULL
  }
  if (stationary) {
    d <- D <- 0
  }
  if (m == 1) {
    D <- max.P <- max.Q <- 0
  }
  else if (is.na(D) & length(xx) <= 2 * m) {
    D <- 0
  }
  else if (is.na(D)) {
    D <- do.call("nsdiffs", c(list(xx, test = seasonal.test,
                                   max.D = max.D), seasonal.test.args))
    if (D > 0 && !is.null(xregg)) {
      diffxreg <- diff(xregg, differences = D, lag = m)
      if (any(apply(diffxreg, 2, is.constant))) {
        D <- D - 1
      }
    }
    if (D > 0) {
      dx <- diff(xx, differences = D, lag = m)
      if (all(is.na(dx)))
        D <- D - 1
    }
  }
  if (D > 0) {
    dx <- diff(xx, differences = D, lag = m)
  }
  else {
    dx <- xx
  }
  if (!is.null(xregg)) {
    if (D > 0) {
      diffxreg <- diff(xregg, differences = D, lag = m)
    }
    else {
      diffxreg <- xregg
    }
  }
  if (is.na(d)) {
    d <- do.call("ndiffs", c(list(dx, test = test, max.d = max.d),
                             test.args))
    if (d > 0 && !is.null(xregg)) {
      diffxreg <- diff(diffxreg, differences = d, lag = 1)
      if (any(apply(diffxreg, 2, is.constant))) {
        d <- d - 1
      }
    }
    if (d > 0) {
      diffdx <- diff(dx, differences = d, lag = 1)
      if (all(is.na(diffdx)))
        d <- d - 1
    }
  }
  if (D >= 2) {
    warning("Having more than one seasonal differences is not recommended. Please consider using only one seasonal difference.")
  }
  else if (D + d > 2) {
    warning("Having 3 or more differencing operations is not recommended. Please consider reducing the total number of differences.")
  }
  if (d > 0) {
    dx <- diff(dx, differences = d, lag = 1)
  }
  if (length(dx) == 0L)
    stop("Not enough data to proceed")
  else if (is.constant(dx)) {
    if (is.null(xreg)) {
      if (D > 0 && d == 0) {
        fit <- arima2(x, order = c(0, d, 0), seasonal = list(order = c(0,
                                                                      D, 0), period = m), include.mean = TRUE,
                     fixed = mean(dx/m, na.rm = TRUE), method = method, nrestart=nrestart,
                     ...)
      }
      else if (D > 0 && d > 0) {
        fit <- arima2(x, order = c(0, d, 0), seasonal = list(order = c(0,
                                                                      D, 0), period = m), method = method,nrestart=nrestart, ...)
      }
      else if (d == 2) {
        fit <- arima2(x, order = c(0, d, 0), method = method, nrestart=nrestart,
                     ...)
      }
      else if (d < 2) {
        fit <- arima2(x, order = c(0, d, 0), include.mean = TRUE,
                     fixed = mean(dx, na.rm = TRUE), method = method, nrestart=nrestart,
                     ...)
      }
      else {
        stop("Data follow a simple polynomial and are not suitable for ARIMA modelling.")
      }
    }
    else {
      if (D > 0) {
        fit <- arima2(x, order = c(0, d, 0), seasonal = list(order = c(0,
                                                                      D, 0), period = m), xreg = xreg, method = method, nrestart=nrestart,
                     ...)
      }
      else {
        fit <- arima2(x, order = c(0, d, 0), xreg = xreg,
                     method = method, nrestart=nrestart, ...)
      }
    }
    fit$x <- orig.x
    fit$series <- series
    fit$call <- match.call()
    fit$call$x <- data.frame(x = x)
    return(fit)
  }
  if (m > 1) {
    if (max.P > 0) {
      max.p <- min(max.p, m - 1)
    }
    if (max.Q > 0) {
      max.q <- min(max.q, m - 1)
    }
  }
  if (approximation) {
    if (!is.null(truncate)) {
      tspx <- stats::tsp(x)
      if (length(x) > truncate) {
        x <- stats::ts(utils::tail(x, truncate), end = tspx[2], frequency = tspx[3])
      }
    }
    if (D == 0) {
      fit <- try(arima2::arima2(x, order = c(0, d, 0), xreg = xreg, nrestart=nrestart,
                              ...), silent = TRUE)
    }
    else {
      fit <- try(arima2::arima2(x, order = c(0, d, 0), seasonal = list(order = c(0,
                                                                               D, 0), period = m), xreg = xreg, nrestart=nrestart, ...), silent = TRUE)
    }
    if (!is.element("try-error", class(fit))) {
      offset <- -2 * fit$loglik - serieslength * log(fit$sigma2)
    }
    else {
      offset <- 0
    }
  }
  else {
    offset <- 0
  }
  allowdrift <- allowdrift & (d + D) == 1
  allowmean <- allowmean & (d + D) == 0
  constant <- allowdrift | allowmean
  ### NEEDS WORK#####################################################################################################
  if (approximation && trace) {
    cat("\n Fitting models using approximations to speed things up...\n")
  }
  if (!stepwise) {
    bestfit <- search.arima(x, d, D, max.p, max.q, max.P,
                            max.Q, max.order, stationary, ic, trace, approximation,
                            method = method, xreg = xreg, offset = offset, allowdrift = allowdrift,
                            allowmean = allowmean, parallel = parallel, num.cores = num.cores,
                            ...)
    bestfit$call <- match.call()
    bestfit$call$x <- data.frame(x = x)
    bestfit$lambda <- lambda
    bestfit$x <- orig.x
    bestfit$series <- series
    bestfit$fitted <- fitted.Arima(bestfit)
    if (trace) {
      cat("\n\n Best model:", arima.string(bestfit, padding = TRUE),
          "\n\n")
    }
    return(bestfit)
  }
  if (length(x) < 10L) {
    start.p <- min(start.p, 1L)
    start.q <- min(start.q, 1L)
    start.P <- 0L
    start.Q <- 0L
  }
  p <- start.p <- min(start.p, max.p)
  q <- start.q <- min(start.q, max.q)
  P <- start.P <- min(start.P, max.P)
  Q <- start.Q <- min(start.Q, max.Q)
  results <- matrix(NA, nrow = nmodels, ncol = 8)
  bestfit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                                                         D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart, ...)
  results[1, ] <- c(p, d, q, P, D, Q, constant, bestfit$ic)
  fit <- myarima(x, order = c(0, d, 0), seasonal = c(0, D,
                                                     0), constant = constant, ic, trace, approximation, method = method,
                 offset = offset, xreg = xreg,nrestart=nrestart, ...)
  results[2, ] <- c(0, d, 0, 0, D, 0, constant, fit$ic)
  if (fit$ic < bestfit$ic) {
    bestfit <- fit
    p <- q <- P <- Q <- 0
  }
  k <- 2
  if (max.p > 0 || max.P > 0) {
    fit <- myarima(x, order = c(max.p > 0, d, 0), seasonal = c((m >
                                                                  1) & (max.P > 0), D, 0), constant = constant, ic,
                   trace, approximation, method = method, offset = offset,
                   xreg = xreg,nrestart=nrestart, ...)
    results[k + 1, ] <- c(max.p > 0, d, 0, (m > 1) & (max.P >
                                                        0), D, 0, constant, fit$ic)
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- (max.p > 0)
      P <- (m > 1) & (max.P > 0)
      q <- Q <- 0
    }
    k <- k + 1
  }
  if (max.q > 0 || max.Q > 0) {
    fit <- myarima(x, order = c(0, d, max.q > 0), seasonal = c(0,
                                                               D, (m > 1) & (max.Q > 0)), constant = constant,
                   ic, trace, approximation, method = method, offset = offset,
                   xreg = xreg,nrestart=nrestart, ...)
    results[k + 1, ] <- c(0, d, max.q > 0, 0, D, (m > 1) &
                            (max.Q > 0), constant, fit$ic)
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- P <- 0
      Q <- (m > 1) & (max.Q > 0)
      q <- (max.q > 0)
    }
    k <- k + 1
  }
  if (constant) {
    fit <- myarima(x, order = c(0, d, 0), seasonal = c(0,
                                                       D, 0), constant = FALSE, ic, trace, approximation,
                   method = method, offset = offset, xreg = xreg,nrestart=nrestart, ...)
    results[k + 1, ] <- c(0, d, 0, 0, D, 0, 0, fit$ic)
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- q <- P <- Q <- 0
    }
    k <- k + 1
  }
  startk <- 0
  while (startk < k && k < nmodels) {
    startk <- k
    if (P > 0 && newmodel(p, d, q, P - 1, D, Q, constant,
                          results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P -
                                                           1, D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p, d, q, P - 1, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        P <- (P - 1)
        next
      }
    }
    if (Q > 0 && newmodel(p, d, q, P, D, Q - 1, constant,
                          results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                                                         D, Q - 1), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p, d, q, P, D, Q - 1, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q - 1)
        next
      }
    }
    if (P < max.P && newmodel(p, d, q, P + 1, D, Q, constant,
                              results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P +
                                                           1, D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p, d, q, P + 1, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        P <- (P + 1)
        next
      }
    }
    if (Q < max.Q && newmodel(p, d, q, P, D, Q + 1, constant,
                              results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                                                         D, Q + 1), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p, d, q, P, D, Q + 1, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q + 1)
        next
      }
    }
    if (Q > 0 && P > 0 && newmodel(p, d, q, P - 1, D, Q -
                                   1, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P -
                                                           1, D, Q - 1), constant = constant, ic, trace,
                     approximation, method = method, offset = offset,
                     xreg = xreg,nrestart=nrestart, ...)
      results[k, ] <- c(p, d, q, P - 1, D, Q - 1, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q - 1)
        P <- (P - 1)
        next
      }
    }
    if (Q < max.Q && P > 0 && newmodel(p, d, q, P - 1, D,
                                       Q + 1, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P -
                                                           1, D, Q + 1), constant = constant, ic, trace,
                     approximation, method = method, offset = offset,
                     xreg = xreg,nrestart=nrestart, ...)
      results[k, ] <- c(p, d, q, P - 1, D, Q + 1, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q + 1)
        P <- (P - 1)
        next
      }
    }
    if (Q > 0 && P < max.P && newmodel(p, d, q, P + 1, D,
                                       Q - 1, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P +
                                                           1, D, Q - 1), constant = constant, ic, trace,
                     approximation, method = method, offset = offset,
                     xreg = xreg,nrestart=nrestart, ...)
      results[k, ] <- c(p, d, q, P + 1, D, Q - 1, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q - 1)
        P <- (P + 1)
        next
      }
    }
    if (Q < max.Q && P < max.P && newmodel(p, d, q, P +
                                           1, D, Q + 1, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q), seasonal = c(P +
                                                           1, D, Q + 1), constant = constant, ic, trace,
                     approximation, method = method, offset = offset,
                     xreg = xreg,nrestart=nrestart, ...)
      results[k, ] <- c(p, d, q, P + 1, D, Q + 1, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        Q <- (Q + 1)
        P <- (P + 1)
        next
      }
    }
    if (p > 0 && newmodel(p - 1, d, q, P, D, Q, constant,
                          results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p - 1, d, q), seasonal = c(P,
                                                             D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p - 1, d, q, P, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- (p - 1)
        next
      }
    }
    if (q > 0 && newmodel(p, d, q - 1, P, D, Q, constant,
                          results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q - 1), seasonal = c(P,
                                                             D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p, d, q - 1, P, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        next
      }
    }
    if (p < max.p && newmodel(p + 1, d, q, P, D, Q, constant,
                              results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p + 1, d, q), seasonal = c(P,
                                                             D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p + 1, d, q, P, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- (p + 1)
        next
      }
    }
    if (q < max.q && newmodel(p, d, q + 1, P, D, Q, constant,
                              results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p, d, q + 1), seasonal = c(P,
                                                             D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p, d, q + 1, P, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        next
      }
    }
    if (q > 0 && p > 0 && newmodel(p - 1, d, q - 1, P, D,
                                   Q, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p - 1, d, q - 1), seasonal = c(P,
                                                                 D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p - 1, d, q - 1, P, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        p <- (p - 1)
        next
      }
    }
    if (q < max.q && p > 0 && newmodel(p - 1, d, q + 1,
                                       P, D, Q, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p - 1, d, q + 1), seasonal = c(P,
                                                                 D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p - 1, d, q + 1, P, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        p <- (p - 1)
        next
      }
    }
    if (q > 0 && p < max.p && newmodel(p + 1, d, q - 1,
                                       P, D, Q, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p + 1, d, q - 1), seasonal = c(P,
                                                                 D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p + 1, d, q - 1, P, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        p <- (p + 1)
        next
      }
    }
    if (q < max.q && p < max.p && newmodel(p + 1, d, q +
                                           1, P, D, Q, constant, results[1:k, ])) {
      k <- k + 1
      if (k > nmodels)
        next
      fit <- myarima(x, order = c(p + 1, d, q + 1), seasonal = c(P,
                                                                 D, Q), constant = constant, ic, trace, approximation,
                     method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                     ...)
      results[k, ] <- c(p + 1, d, q + 1, P, D, Q, constant,
                        fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        p <- (p + 1)
        next
      }
    }
    if (allowdrift || allowmean) {
      if (newmodel(p, d, q, P, D, Q, !constant, results[1:k,
      ])) {
        k <- k + 1
        if (k > nmodels)
          next
        fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                                                           D, Q), constant = !constant, ic, trace, approximation,
                       method = method, offset = offset, xreg = xreg,nrestart=nrestart,
                       ...)
        results[k, ] <- c(p, d, q, P, D, Q, !constant,
                          fit$ic)
        if (fit$ic < bestfit$ic) {
          bestfit <- fit
          constant <- !constant
        }
      }
    }
  }
  if (k > nmodels) {
    warning(sprintf("Stepwise search was stopped early due to reaching the model number limit: `nmodels = %i`",
                    nmodels))
  }
  if (approximation && !is.null(bestfit$arma)) {
    if (trace) {
      cat("\n\n Now re-fitting the best model(s) without approximations...\n")
    }
    icorder <- order(results[, 8])
    nmodels <- sum(!is.na(results[, 8]))
    for (i in seq(nmodels)) {
      k <- icorder[i]
      fit <- myarima(x, order = c(results[k, 1], d, results[k,
                                                            3]), seasonal = c(results[k, 4], D, results[k,
                                                                                                        6]), constant = results[k, 7] == 1, ic, trace,
                     approximation = FALSE, method = method, xreg = xreg,nrestart=nrestart,
                     ...)
      if (fit$ic < Inf) {
        bestfit <- fit
        break
      }
    }
  }
  if (bestfit$ic == Inf && !isTRUE(method == "CSS")) {
    if (trace) {
      cat("\n")
    }
    stop("No suitable ARIMA model found")
  }
  bestfit$x <- orig.x
  bestfit$series <- series
  bestfit$ic <- NULL
  bestfit$call <- match.call()
  bestfit$call$x <- data.frame(x = x)
  bestfit$lambda <- lambda
  bestfit$fitted <- fitted.Arima(bestfit)
  if (trace) {
    cat("\n\n Best model:", arima.string(bestfit, padding = TRUE),
        "\n\n")
  }
  return(bestfit)
}
# Calls arima2 from arima2 package and adds data to the returned object
# Also allows refitting to new data
# and drift terms to be included.
myarima <- function(x, order = c(0, 0, 0), seasonal = c(0, 0, 0), constant=TRUE, ic="aic", trace=FALSE, approximation=FALSE, offset=0, xreg=NULL, method = NULL, nrestart=10, ...) {
  # Length of non-missing interior
  missing <- is.na(x)
  firstnonmiss <- utils::head(which(!missing),1)
  lastnonmiss <- utils::tail(which(!missing),1)
  n <- sum(!missing[firstnonmiss:lastnonmiss])
  m <- stats::frequency(x)
  use.season <- (sum(seasonal) > 0) & m > 0
  diffs <- order[2] + seasonal[2]
  if(is.null(method)){
    if (approximation) {
      method <- "CSS"
    } else {
      method <- "CSS-ML"
    }
  }
  if (diffs == 1 && constant) {
    xreg <- `colnames<-`(cbind(drift = 1:length(x), xreg),
                         make.unique(c("drift", if(is.null(colnames(xreg)) && !is.null(xreg)) rep("", NCOL(xreg)) else colnames(xreg))))
    if (use.season) {
      suppressWarnings(fit <- try(arima2::arima2(x = x, order = order, seasonal = list(order = seasonal, period = m), xreg = xreg, method = method, nrestart=nrestart, ...), silent = TRUE))
    } else {
      suppressWarnings(fit <- try(arima2::arima2(x = x, order = order, xreg = xreg, method = method, nrestart=nrestart, ...), silent = TRUE))
    }
  }
  else {
    if (use.season) {
      suppressWarnings(fit <- try(arima2::arima2(x = x, order = order, seasonal = list(order = seasonal, period = m), include.mean = constant, method = method, xreg = xreg, nrestart=nrestart, ...), silent = TRUE))
    } else {
      suppressWarnings(fit <- try(arima2::arima2(x = x, order = order, include.mean = constant, method = method, xreg = xreg, nrestart=nrestart, ...), silent = TRUE))
    }
  }
  if (is.null(xreg)) {
    nxreg <- 0
  } else {
    nxreg <- ncol(as.matrix(xreg))
  }
  if (!is.element("try-error", class(fit))) {
    nstar <- n - order[2] - seasonal[2] * m
    if (diffs == 1 && constant) {
      # fitnames <- names(fit$coef)
      # fitnames[length(fitnames)-nxreg] <- "drift"
      # names(fit$coef) <- fitnames
      fit$xreg <- xreg
    }
    npar <- length(fit$coef[fit$mask]) + 1
    if (method == "CSS") {
      fit$aic <- offset + nstar * log(fit$sigma2) + 2 * npar
    }
    if (!is.na(fit$aic)) {
      fit$bic <- fit$aic + npar * (log(nstar) - 2)
      fit$aicc <- fit$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
      fit$ic <- switch(ic, bic = fit$bic, aic = fit$aic, aicc = fit$aicc)
    }
    else {
      fit$aic <- fit$bic <- fit$aicc <- fit$ic <- Inf
    }
    # Adjust residual variance to be unbiased
    fit$sigma2 <- sum(fit$residuals ^ 2, na.rm = TRUE) / (nstar - npar + 1)

    # Check for unit roots
    minroot <- 2
    if (order[1] + seasonal[1] > 0) {
      testvec <- fit$model$phi
      k <- abs(testvec) > 1e-8
      if (sum(k) > 0) {
        last.nonzero <- max(which(k))
      } else {
        last.nonzero <- 0
      }
      if (last.nonzero > 0) {
        testvec <- testvec[1:last.nonzero]
        proots <- try(polyroot(c(1,-testvec)))
        if (!is.element("try-error", class(proots))) {
          minroot <- min(minroot, abs(proots))
        }
        else fit$ic <- Inf
      }
    }
    if (order[3] + seasonal[3] > 0 & fit$ic < Inf) {
      testvec <- fit$model$theta
      k <- abs(testvec) > 1e-8
      if (sum(k) > 0) {
        last.nonzero <- max(which(k))
      } else {
        last.nonzero <- 0
      }
      if (last.nonzero > 0) {
        testvec <- testvec[1:last.nonzero]
        proots <- try(polyroot(c(1,testvec)))
        if (!is.element("try-error", class(proots))) {
          minroot <- min(minroot, abs(proots))
        }
        else fit$ic <- Inf
      }
    }
    # Avoid bad models
    if (minroot < 1 + 1e-2 | checkarima(fit)) {
      fit$ic <- Inf
    }

    fit$xreg <- xreg

    if (trace) {
      cat("\n", arima.string(fit, padding = TRUE), ":", fit$ic)
    }

    return(structure(fit, class = c("forecast_ARIMA", "ARIMA", "Arima")))
  }
  else {
    # Catch errors due to unused arguments
    if (length(grep("unused argument", fit)) > 0L) {
      stop(fit[1])
    }

    if (trace) {
      cat("\n ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
      if (use.season) {
        cat("(", seasonal[1], ",", seasonal[2], ",", seasonal[3], ")[", m, "]", sep = "")
      }
      if (constant && (order[2] + seasonal[2] == 0)) {
        cat(" with non-zero mean")
      } else if (constant && (order[2] + seasonal[2] == 1)) {
        cat(" with drift        ")
      } else if (!constant && (order[2] + seasonal[2] == 0)) {
        cat(" with zero mean    ")
      } else {
        cat("                   ")
      }
      cat(" :", Inf)
    }
    return(list(ic = Inf))
  }
}

newmodel <- function(p, d, q, P, D, Q, constant, results) {
  n <- nrow(results)
  for (i in 1:n) {
    if(!all(is.na(results[i, seq(7)]))) {
      if (all(c(p, d, q, P, D, Q, constant) == results[i, 1:7])) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

arima.string <- function(object, padding=FALSE) {
  order <- object$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
  if (m > 1 && sum(order[4:6]) > 0) {
    result <- paste(result, "(", order[4], ",", order[5], ",", order[6], ")[", m, "]", sep = "")
  }
  if (padding && m > 1 && sum(order[4:6]) == 0) {
    result <- paste(result, "         ", sep = "")
    if (m <= 9) {
      result <- paste(result, " ", sep = "")
    } else if (m <= 99) {
      result <- paste(result, "  ", sep = "")
    } else {
      result <- paste(result, "   ", sep = "")
    }
  }
  if (!is.null(object$xreg)) {
    if (NCOL(object$xreg) == 1 && is.element("drift", names(object$coef))) {
      result <- paste(result, "with drift        ")
    } else {
      result <- paste("Regression with", result, "errors")
    }
  }
  else {
    if (is.element("constant", names(object$coef)) || is.element("intercept", names(object$coef))) {
      result <- paste(result, "with non-zero mean")
    } else if (order[2] == 0 && order[5] == 0) {
      result <- paste(result, "with zero mean    ")
    } else {
      result <- paste(result, "                  ")
    }
  }
  if (!padding) {
    # Strip trailing spaces
    result <- gsub("[ ]*$", "", result)
  }
  return(result)
}

#' @export
summary.Arima <- function(object, ...) {
  class(object) <- c("summary.Arima", class(object))
  object
}

# Check that Arima object has positive coefficient variances without returning warnings
checkarima <- function(object) {
  suppressWarnings(test <- any(is.nan(sqrt(diag(object$var.coef)))))
  return(test)
}

#' Is an object constant?
#'
#' Returns true if the object's numerical values do not vary.
#'
#'
#' @param x object to be tested
#' @export
is.constant <- function(x) {
  x <- as.numeric(x)
  y <- rep(x[1], length(x))
  return(isTRUE(all.equal(x, y)))
}

search.arima <- function(x, d=NA, D=NA, max.p=5, max.q=5,
                         max.P=2, max.Q=2, max.order=5, stationary=FALSE, ic=c("aic", "aicc", "bic"),
                         trace=FALSE, approximation=FALSE, xreg=NULL, offset=offset, allowdrift=TRUE,
                         allowmean=TRUE, parallel=FALSE, num.cores=2, ...) {
  # dataname <- substitute(x)
  ic <- match.arg(ic)
  m <- stats::frequency(x)

  allowdrift <- allowdrift & (d + D) == 1
  allowmean <- allowmean & (d + D) == 0

  maxK <- (allowdrift | allowmean)

  # Choose model orders
  # Serial - technically could be combined with the code below
  if (parallel == FALSE) {
    best.ic <- Inf
    for (i in 0:max.p) {
      for (j in 0:max.q) {
        for (I in 0:max.P) {
          for (J in 0:max.Q) {
            if (i + j + I + J <= max.order) {
              for (K in 0:maxK) {
                fit <- myarima(
                  x, order = c(i, d, j), seasonal = c(I, D, J),
                  constant = (K == 1), trace = trace, ic = ic, approximation = approximation,
                  offset = offset, xreg = xreg, ...
                )
                if (fit$ic < best.ic) {
                  best.ic <- fit$ic
                  bestfit <- fit
                  constant <- (K == 1)
                }
              }
            }
          }
        }
      }
    }
  } else if (parallel == TRUE) {
    to.check <- WhichModels(max.p, max.q, max.P, max.Q, maxK)

    par.all.arima <- function(l, max.order) {
      .tmp <- UndoWhichModels(l)
      i <- .tmp[1]
      j <- .tmp[2]
      I <- .tmp[3]
      J <- .tmp[4]
      K <- .tmp[5] == 1
      if (i + j + I + J <= max.order) {
        fit <- myarima(
          x, order = c(i, d, j), seasonal = c(I, D, J), constant = (K == 1),
          trace = trace, ic = ic, approximation = approximation, offset = offset, xreg = xreg,
          ...
        )
      }
      if (exists("fit")) {
        return(cbind(fit, K))
      } else {
        return(NULL)
      }
    }

    if (is.null(num.cores)) {
      num.cores <- parallel::detectCores()
    }

    all.models <- parallel::mclapply(X = to.check, FUN = par.all.arima, max.order=max.order, mc.cores = num.cores)

    # Removing null elements
    all.models <- all.models[!sapply(all.models, is.null)]

    # Choosing best model
    best.ic <- Inf
    for (i in 1:length(all.models)) {
      if (!is.null(all.models[[i]][, 1]$ic) && all.models[[i]][, 1]$ic < best.ic) {
        bestfit <- all.models[[i]][, 1]
        best.ic <- bestfit$ic
        constant <- unlist(all.models[[i]][1, 2])
      }
    }
    class(bestfit) <- c("forecast_ARIMA", "ARIMA", "Arima")
  }

  if (exists("bestfit")) {
    # Refit using ML if approximation used for IC
    if (approximation) {
      if (trace) {
        cat("\n\n Now re-fitting the best model(s) without approximations...\n")
      }
      # constant <- length(bestfit$coef) - ncol(xreg) > sum(bestfit$arma[1:4])
      newbestfit <- myarima(
        x, order = bestfit$arma[c(1, 6, 2)],
        seasonal = bestfit$arma[c(3, 7, 4)], constant = constant, ic,
        trace = FALSE, approximation = FALSE, xreg = xreg, ...
      )
      if (newbestfit$ic == Inf) {
        # Final model is lousy. Better try again without approximation
        # warning("Unable to fit final model using maximum likelihood. AIC value approximated")
        bestfit <- search.arima(
          x, d = d, D = D, max.p = max.p, max.q = max.q,
          max.P = max.P, max.Q = max.Q, max.order = max.order, stationary = stationary,
          ic = ic, trace = trace, approximation = FALSE, xreg = xreg, offset = offset,
          allowdrift = allowdrift, allowmean = allowmean,
          parallel = parallel, num.cores = num.cores, ...
        )
        bestfit$ic <- switch(ic, bic = bestfit$bic, aic = bestfit$aic, aicc = bestfit$aicc)
      }
      else {
        bestfit <- newbestfit
      }
    }
  }
  else {
    stop("No ARIMA model able to be estimated")
  }

  bestfit$x <- x
  bestfit$series <- deparse(substitute(x))
  bestfit$ic <- NULL
  bestfit$call <- match.call()

  if (trace) {
    cat("\n\n")
  }

  return(bestfit)
}

# Set up seasonal dummies using Fourier series
SeasDummy <- function(x) {
  n <- length(x)
  m <- stats::frequency(x)
  if (m == 1) {
    stop("Non-seasonal data")
  }
  tt <- 1:n
  fmat <- matrix(NA, nrow = n, ncol = 2 * m)
  for (i in 1:m) {
    fmat[, 2 * i] <- sin(2 * pi * i * tt / m)
    fmat[, 2 * (i - 1) + 1] <- cos(2 * pi * i * tt / m)
  }
  return(fmat[, 1:(m - 1)])
}

# CANOVA-HANSEN TEST
# Largely based on uroot package code for CH.test()
SD.test <- function(wts, s=stats::frequency(wts)) {
  if (any(is.na(wts))) {
    stop("Series contains missing values. Please choose order of seasonal differencing manually.")
  }
  if (s == 1) {
    stop("Not seasonal data")
  }
  t0 <- stats::start(wts)
  N <- length(wts)
  if (N <= s) {
    stop("Insufficient data")
  }
  frec <- rep(1, as.integer((s + 1) / 2))
  ltrunc <- round(s * (N / 100) ^ 0.25)
  R1 <- as.matrix(SeasDummy(wts))
  lmch <- stats::lm(wts ~ R1, na.action = stats::na.exclude) # run the regression : y(i)=mu+f(i)'gamma(i)+e(i)
  Fhat <- Fhataux <- matrix(nrow = N, ncol = s - 1)
  for (i in 1:(s - 1))
    Fhataux[, i] <- R1[, i] * stats::residuals(lmch)
  for (i in 1:N) {
    for (n in 1:(s - 1))
      Fhat[i, n] <- sum(Fhataux[1:i, n])
  }
  wnw <- 1 - seq(1, ltrunc, 1) / (ltrunc + 1)
  Ne <- nrow(Fhataux)
  Omnw <- 0
  for (k in 1:ltrunc)
    Omnw <- Omnw + (t(Fhataux)[, (k + 1):Ne] %*% Fhataux[1:(Ne - k), ]) * wnw[k]
  Omfhat <- (crossprod(Fhataux) + Omnw + t(Omnw)) / Ne
  sq <- seq(1, s - 1, 2)
  frecob <- rep(0, s - 1)
  for (i in 1:length(frec)) {
    if (frec[i] == 1 && i == as.integer(s / 2)) {
      frecob[sq[i]] <- 1
    }
    if (frec[i] == 1 && i < as.integer(s / 2)) {
      frecob[sq[i]] <- frecob[sq[i] + 1] <- 1
    }
  }
  a <- length(which(frecob == 1))
  A <- matrix(0, nrow = s - 1, ncol = a)
  j <- 1
  for (i in 1:(s - 1)) {
    if (frecob[i] == 1) {
      A[i, j] <- 1
      ifelse(frecob[i] == 1, j <- j + 1, j <- j)
    }
  }
  tmp <- t(A) %*% Omfhat %*% A
  problems <- (min(svd(tmp)$d) < .Machine$double.eps)
  if (problems) {
    stL <- 0
  } else {
    stL <- (1 / N ^ 2) * sum(diag(solve(tmp, tol = 1e-25) %*% t(A) %*% t(Fhat) %*% Fhat %*% A))
  }
  return(stL)
}
