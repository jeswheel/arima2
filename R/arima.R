#' ARIMA Modeling of Time Series
#'
#' Fit an ARIMA model to a univariate time series. This function builds on
#' the ARIMA model fitting approach used in [stats::arima()] by fitting
#' model parameters via a random restart algorithm.
#'
#' @param diffuseControl Boolean indicator of whether or initial observations
#'    will have likelihood values ignored if controlled by the diffuse prior,
#'    i.e., have a Kalman gain of at least 1e4.
#' @param max_iters Maximum number of random restarts for methods "CSS-ML" and
#'    "ML". If set to 1, the results of this algorithm is the same as
#'    [stats::arima()] if argument \code{diffuseControl} is also set as TRUE.
#'    \code{max_iters} is often not reached because the condition
#'    \code{max_repeats} is typically achieved first.
#' @param max_repeats Integer. If the last \code{max_repeats} random starts did
#'    not result in improved likelihoods, then stop the search. Each result of
#'    the optim function is only considered to improve the likelihood if it does
#'    so by more than \code{eps_tol}.
#' @param max_inv_root positive numeric value less than or equal to 1. This
#'    number represents the maximum size of the inverted
#'    MA or AR polynomial roots for a new parameter estimate to be considered an
#'    improvement to previous estimates. Concerns of numeric stability arise
#'    when the size of polynomial roots are near unity circle. The default value
#'    1 means that the the parameter values corresponding with the best
#'    log-likelihood will be returned, even if they are near unity.
#'    Suitable values of this parameter are near the value 1.
#' @param min_inv_root_dist positive numeric value less than 1. This number
#'    represents the minimum distance between AR and MA polynomial roots for a
#'    new parameter estimate to be considered an improvement on previous
#'    estimates. This is intended to avoid the possibility of returning
#'    parameter estimates with nearly canceling roots. Appropriate choices are
#'    values near 0.
#' @param eps_tol Tolerance for accepting a new solution to be better than a
#'    previous solution in terms of log-likelihood. The default corresponds to a
#'    one ten-thousandth unit increase in log-likelihood.
#' @param init_method Method used to randomly sample parameter initializations.
#'    \code{init_method = "DL"} will sample parameters using the Durbin-Levinson
#'    algorithm, described by Monahan (1984). If
#'    \code{init_method = "UnifRoots"}, then inverted roots of AR and MA
#'    polynomials will be sampled uniformly from the complex unit circle.
#' @inheritParams stats::arima
#' @returns
#' A list of class \code{c("Arima2", "Arima")}. This list contains all of the
#' same elements as the output of [stats::arima], along with some additional
#' elements.  All elements of the output list are:
#' \describe{
#'    \item{`coef`}{A vector of AR, MA, and regression coefficients. These can
#'    be extracted by the [stats::coef] method.}
#'    \item{`sigma2`}{The MLE of the variance of the innovations.}
#'    \item{`var.coef`}{The estimated variance matrix of the coefficients
#'     `coef`, which can be extracted by the [stats::vcov] method.}
#'    \item{`mask`}{A vector containing boolean values, indicating which
#'     parameters of the model were estimated.}
#'    \item{`loglik`}{The maximized log-likelihood (of the differenced data).}
#'    \item{`aic`}{The AIC value corresponding to the log-likelihood.}
#'    \item{`arma`}{A compact form of the model specification, as a vector
#'     giving the number of AR, MA, seasonal AR and seasonal MA coefficients,
#'     plus the period and the number of non-seasonal and seasonal differences.}
#'    \item{`residuals`}{The fitted innovations.}
#'    \item{`call`}{The matched call.}
#'    \item{`series`}{The name of the series x.}
#'    \item{`code`}{The convergence value returned by [stats::optim].}
#'    \item{`n.cond`}{The number of initial observations not used in the
#'     fitting.}
#'    \item{`nobs`}{The number of observations used for the fitting.}
#'    \item{`model`}{A list representing the Kalman Filter used in the fitting.}
#'    \item{`x`}{The input time series.}
#'    \item{`num_starts`}{Number of restarts before convergence criteria was
#'    satisfied.}
#'    \item{`all_values`}{Numeric vector of length `num_starts` containing the
#'    loglikelihood of every parameter initialization.}
#' }
#'
#' @references Monahan, John F. (1984) A note on enforcing stationarity in autoregressive-moving average models. \emph{Biometrika}, \bold{71}(2), 403--404.
#'
#' @export
#' @examples
#' # example code
#' set.seed(12345)
#' arima(miHuron_level$Average, order = c(2, 0, 1), max_iters = 100)
#'
#' @useDynLib arima2, ARIMA_transPars, ARIMA_CSS, TSconv, getQ0, getQ0bis, ARIMA_Like, ARIMA_Invtrans, ARIMA_Gradtrans, ARIMA_undoPars, .fixes = "C_"
arima <- function(x, order = c(0L, 0L, 0L),
                  seasonal = list(order = c(0L, 0L, 0L), period = NA),
                  xreg = NULL, include.mean = TRUE,
                  transform.pars = TRUE, fixed = NULL, init = NULL,
                  method = c("CSS-ML", "ML", "CSS"),
                  init_method = c("DL", "UnifRoots"),
                  n.cond,
                  SSinit = c("Rossignol2011", "Gardner1980"),
                  optim.method = "BFGS",
                  optim.control = list(), kappa = 1e6,
                  diffuseControl = TRUE,
                  max_iters = 100,
                  max_repeats = 10,
                  max_inv_root = 1,
                  min_inv_root_dist = 0,
                  eps_tol = 1e-4)
{

  #  This function is based on the arima function of the stats package
  #  of R. Below the copright statement of the arima function is reproduced.
  #
  #       File src/library/stats/R/arima.R
  #       Part of the stats package
  #
  #  Copyright (C) 2002-16 The R Core Team
  #
  #  This program is free software; you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation; either version 2 of the License, or
  #  (at your option) any later version.
  #
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.
  #
  #  A copy of the GNU General Public License is available at
  #  http://www.r-project.org/Licenses/

  #
  # (arima2) Date: Dec 6, 2022
  # Revised: Jan 23, 2023

  "%+%" <- function(a, b) .Call(C_TSconv, a, b)

  SSinit <- match.arg(SSinit)
  init_method <- match.arg(init_method)
  SS.G <- SSinit == "Gardner1980"
  ## helper of armafn(), called by optim()
  upARIMA <- function(mod, phi, theta) {
    p <- length(phi); q <- length(theta)
    mod$phi <- phi; mod$theta <- theta
    r <- max(p, q + 1L)
    if(p > 0) mod$T[1L:p, 1L] <- phi
    if(r > 1L)
      mod$Pn[1L:r, 1L:r] <-
      if(SS.G) .Call(C_getQ0, phi, theta)
    else .Call(C_getQ0bis, phi, theta, tol = 0)# tol=0: less checking
    else
      mod$Pn[1L, 1L] <- if (p > 0) 1/(1 - phi^2) else 1
    mod$a[] <- 0
    mod
  }

  arimaSS <- function(y, mod) {
    ## next call changes mod components a, P, Pn so beware!
    .Call(C_ARIMA_Like, y, mod, 0L, TRUE, diffuseControl)
  }

  ## the objective function called by optim()
  armafn <- function(p, trans) {
    par <- coef
    par[mask] <- p
    trarma <- .Call(C_ARIMA_transPars, par, arma, trans)
    if(is.null(Z <- tryCatch(upARIMA(mod, trarma[[1L]], trarma[[2L]]),
                             error = function(e) NULL)))
      return(.Machine$double.xmax)# bad parameters giving error, e.g. in solve(.)
    if(ncxreg > 0) x <- x - xreg %*% par[narma + (1L:ncxreg)]
    ## next call changes Z components a, P, Pn so beware!
    res <- .Call(C_ARIMA_Like, x, Z, 0L, FALSE, diffuseControl)
    s2 <- res[1L]/res[3L]
    0.5*(log(s2) + res[2L]/res[3L])
  }

  armaCSS <- function(p) {
    par <- as.double(fixed)
    par[mask] <- p
    trarma <- .Call(C_ARIMA_transPars, par, arma, FALSE)
    if(ncxreg > 0) x <- x - xreg %*% par[narma + (1L:ncxreg)]
    res <- .Call(C_ARIMA_CSS, x, arma, trarma[[1L]], trarma[[2L]],
                 as.integer(ncond), FALSE)
    0.5 * log(res)
  }

  arCheck <- function(ar) {
    p <- max(which(c(1, -ar) != 0)) - 1
    if(!p) return(TRUE)
    all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
  }

  if (max_inv_root < 0) stop("max_inv_root must be positive.")
  if (min_inv_root_dist > 1 || min_inv_root_dist < 0) stop("min_inv_root_dist must in the interval [0, 1].")

  if (max_inv_root == 1) {
    do_maxroot_test <- FALSE
  } else {
    do_maxroot_test <- TRUE
  }

  if (min_inv_root_dist == 0) {
    do_min_dist_test <- FALSE
  } else {
    do_min_dist_test <- TRUE
  }

  series <- deparse1(substitute(x))
  if(NCOL(x) > 1L)
    stop("only implemented for univariate time series")
  method <- match.arg(method)

  x <- stats::as.ts(x)
  if(!is.numeric(x))
    stop("'x' must be numeric")
  storage.mode(x) <- "double"  # a precaution
  dim(x) <- NULL
  n <- length(x)

  if (!is.logical(diffuseControl)) {
    stop("'diffuseControl' must be logical.")
  }

  if(!missing(order))
    if(!is.numeric(order) || length(order) != 3L || any(order < 0))
      stop("'order' must be a non-negative numeric vector of length 3")
  if(!missing(seasonal))
    if(is.list(seasonal)) {
      if(is.null(seasonal$order))
        stop("'seasonal' must be a list with component 'order'")
      if(!is.numeric(seasonal$order) || length(seasonal$order) != 3L
         || any(seasonal$order < 0L))
        stop("'seasonal$order' must be a non-negative numeric vector of length 3")
    } else if(is.numeric(order)) {
      if(length(order) == 3L) seasonal <- list(order=seasonal)
      else ("'seasonal' is of the wrong length")
    } else stop("'seasonal' must be a list with component 'order'")

  if (is.null(seasonal$period) || is.na(seasonal$period) || seasonal$period == 0)
    seasonal$period <- stats::frequency(x)
  arma <- as.integer(c(order[-2L], seasonal$order[-2L], seasonal$period,
                       order[2L], seasonal$order[2L]))
  narma <- sum(arma[1L:4L])

  xtsp <- stats::tsp(x)
  stats::tsp(x) <- NULL
  Delta <- 1.
  for(i in seq_len(order[2L])) Delta <- Delta %+% c(1., -1.)
  for(i in seq_len(seasonal$order[2L]))
    Delta <- Delta %+% c(1, rep.int(0, seasonal$period-1), -1)
  Delta <- - Delta[-1L]
  nd <- order[2L] + seasonal$order[2L]
  n.used <- sum(!is.na(x)) - length(Delta)
  if (is.null(xreg)) {
    ncxreg <- 0L
  } else {
    nmxreg <- deparse1(substitute(xreg))
    if (NROW(xreg) != n) stop("lengths of 'x' and 'xreg' do not match")
    ncxreg <- NCOL(xreg)
    xreg <- as.matrix(xreg)
    storage.mode(xreg) <- "double"
  }
  class(xreg) <- NULL
  if (ncxreg > 0L && is.null(colnames(xreg)))
    colnames(xreg) <-
    if(ncxreg == 1L) nmxreg else paste0(nmxreg, 1L:ncxreg)
  if (include.mean && (nd == 0L)) {
    xreg <- cbind(intercept = rep(1, n), xreg = xreg)
    ncxreg <- ncxreg + 1L
  }
  if(method == "CSS-ML") {
    anyna <- anyNA(x)
    if(ncxreg) anyna <- anyna || anyNA(xreg)
    if(anyna) method <- "ML"
  }

  if (method == "CSS" || method == "CSS-ML") {
    ncond <- order[2L] + seasonal$order[2L] * seasonal$period
    ncond1 <- order[1L] + seasonal$period * seasonal$order[1L]
    ncond <- ncond + if(!missing(n.cond)) max(n.cond, ncond1) else ncond1
  } else ncond <- 0

  if (is.null(fixed)) fixed <- rep(NA_real_, narma + ncxreg)
  else if(length(fixed) != narma + ncxreg) stop("wrong length for 'fixed'")
  mask <- is.na(fixed)

  # If both 'fixed' and 'init' provided, check that they are
  # the same (this doesn't affect model fitting, but becomes
  # important when checking the stationarity of the fitted model,
  # which only checks the vector init.)
  if (any(mask)) {
    if (is.null(init)) {
      init <- rep(NA_real_, length(fixed))
      init[!mask] <- fixed[!mask]
    } else if (any(fixed[!mask] != init[!mask])) {
      init[!mask] = fixed[!mask]
      warning("Both arguments 'fixed' and 'init' provided, but provided coefficients did not match: setting non-missing 'init' values to corresponding values of 'fixed'.")
    }
  }

  ##    if(!any(mask)) stop("all parameters were fixed")
  no.optim <- !any(mask)
  if(no.optim) transform.pars <- FALSE
  if(transform.pars) {
    ind <- arma[1L] + arma[2L] + seq_len(arma[3L])
    if (any(!mask[seq_len(arma[1L])]) || any(!mask[ind])) {
      warning("some AR parameters were fixed: setting transform.pars = FALSE")
      transform.pars <- FALSE
    }
  }
  init0 <- rep.int(0, narma)
  parscale <- rep(1, narma)
  if (ncxreg) {
    cn <- colnames(xreg)
    orig.xreg <- (ncxreg == 1L) || any(!mask[narma + 1L:ncxreg])
    if (!orig.xreg) {
      S <- svd(stats::na.omit(xreg))
      xreg <- xreg %*% S$v
    }
    dx <- x
    dxreg <- xreg
    if(order[2L] > 0L) {
      dx <- diff(dx, 1L, order[2L])
      dxreg <- diff(dxreg, 1L, order[2L])
    }
    if(seasonal$period > 1L && seasonal$order[2L] > 0) {
      dx <- diff(dx, seasonal$period, seasonal$order[2L])
      dxreg <- diff(dxreg, seasonal$period, seasonal$order[2L])
    }
    fit <- if(length(dx) > ncol(dxreg))
      stats::lm(dx ~ dxreg - 1, na.action = stats::na.omit)
    else list(rank = 0L)
    if(fit$rank == 0L) {
      ## Degenerate model. Proceed anyway so as not to break old code
      fit <- stats::lm(x ~ xreg - 1, na.action = stats::na.omit)
    }
    isna <- is.na(x) | apply(xreg, 1L, anyNA)
    n.used <- sum(!isna) - length(Delta)
    init0 <- c(init0, coef(fit))
    ses <- summary(fit)$coefficients[, 2L]
    parscale <- c(parscale, 10 * ses)
  }
  if (n.used <= 0) stop("too few non-missing observations")

  if(!is.null(init)) {
    if(length(init) != length(init0))
      stop("'init' is of the wrong length")
    if(any(ind <- is.na(init))) init[ind] <- init0[ind]
    if(method == "ML") {
      ## check stationarity
      if(arma[1L] > 0)
        if(!arCheck(init[1L:arma[1L]]))
          stop("non-stationary AR part")
      if(arma[3L] > 0)
        if(!arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
          stop("non-stationary seasonal AR part")
      if(transform.pars)
        init <- .Call(C_ARIMA_Invtrans, as.double(init), arma)
    }
  } else init <- init0

  coef <- as.double(fixed)
  if(!("parscale" %in% names(optim.control)))
    optim.control$parscale <- parscale[mask]

  if(method == "CSS") {
    i_start = NULL
    res <- if(no.optim)
      list(convergence=0L, par=numeric(), value=armaCSS(numeric()))
    else
      stats::optim(init[mask], armaCSS,  method = optim.method, hessian = TRUE,
            control = optim.control)
    if(res$convergence > 0)
      warning(gettextf("possible convergence problem: optim gave code = %d",
                       res$convergence), domain = NA)
    coef[mask] <- res$par
    ## set model for predictions
    trarma <- .Call(C_ARIMA_transPars, coef, arma, FALSE)
    mod <- stats::makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa, SSinit)
    if(ncxreg > 0) x <- x - xreg %*% coef[narma + (1L:ncxreg)]
    arimaSS(x, mod)
    val <- .Call(C_ARIMA_CSS, x, arma, trarma[[1L]], trarma[[2L]],
                 as.integer(ncond), TRUE)
    sigma2 <- val[[1L]]
    var <- if(no.optim) numeric() else solve(res$hessian * n.used)

    value <- 2 * n.used * res$value + n.used + n.used * log(2 * pi)
    all_values <- value
  } else {
    if (method == "CSS-ML") {
      res <- if(no.optim)
        list(convergence=0L, par=numeric(), value=armaCSS(numeric()))
      else
        stats::optim(init[mask], armaCSS,  method = optim.method,
              hessian = FALSE, control = optim.control)
      if(res$convergence == 0) init[mask] <- res$par
      ## check stationarity
      if(arma[1L] > 0)
        if(!arCheck(init[1L:arma[1L]]))
          stop("non-stationary AR part from CSS")
      if(arma[3L] > 0)
        if(!arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
          stop("non-stationary seasonal AR part from CSS")
      ncond <- 0L
    }

    # Start of ML estimation

    best_value <- Inf
    best_coef <- coef
    converged <- FALSE
    all_values <- c()
    num_repeat <- 1

    for (i_start in 1:max_iters) {  # Do random restarts.

      if (i_start == 1) {  # use CSS start as baseline, same result as stats::arima
        if(transform.pars) {
          new_init <- .Call(C_ARIMA_Invtrans, init, arma)
          ## enforce invertibility
          if(arma[2L] > 0) {
            ind <- arma[1L] + 1L:arma[2L]
            new_init[ind] <- .maInvert(new_init[ind])
          }
          if(arma[4L] > 0) {
            ind <- sum(arma[1L:3L]) + 1L:arma[4L]
            new_init[ind] <- .maInvert(new_init[ind])
          }
        } else {
          new_init <- init
        }
        trarma <- .Call(C_ARIMA_transPars, new_init, arma, transform.pars)
        mod <- stats::makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa, SSinit)
        best_res <- if(no.optim)
          list(convergence = 0, par = numeric(),
               value = armafn(numeric(), as.logical(transform.pars)))
        else
          stats::optim(new_init[mask], armafn, method = optim.method,
                       hessian = TRUE, control = optim.control,
                       trans = as.logical(transform.pars))
        if(best_res$convergence > 0)
          warning(gettextf("possible convergence problem: optim gave code = %d",
                           best_res$convergence), domain = NA)
        best_coef[mask] <- best_res$par
        if(transform.pars) {
          ## enforce invertibility
          if(arma[2L] > 0L) {
            ind <- arma[1L] + 1L:arma[2L]
            if(all(mask[ind]))
              best_coef[ind] <- .maInvert(best_coef[ind])
          }
          if(arma[4L] > 0L) {
            ind <- sum(arma[1L:3L]) + 1L:arma[4L]
            if(all(mask[ind]))
              best_coef[ind] <- .maInvert(best_coef[ind])
          }
          if(any(best_coef[mask] != best_res$par))  {  # need to re-fit
            oldcode <- best_res$convergence
            best_res <- stats::optim(best_coef[mask], armafn, method = optim.method,
                                hessian = TRUE,
                                control = list(maxit = 0L,
                                               parscale = optim.control$parscale),
                                trans = TRUE)
            best_res$convergence <- oldcode
            best_coef[mask] <- best_res$par
          }
          ## do it this way to ensure hessian was computed inside
          ## stationarity region
          A <- .Call(C_ARIMA_Gradtrans, as.double(best_coef), arma)
          A <- A[mask, mask]
          best_var <- crossprod(A, solve(best_res$hessian * n.used, A))
          best_coef <- .Call(C_ARIMA_undoPars, best_coef, arma)
        } else best_var <- if(no.optim) numeric() else solve(best_res$hessian * n.used)

        best_trarma <- .Call(C_ARIMA_transPars, best_coef, arma, FALSE)
        best_mod <- stats::makeARIMA(best_trarma[[1L]], best_trarma[[2L]], Delta, kappa, SSinit)
        best_val <- if(ncxreg > 0L)
          arimaSS(x - xreg %*% best_coef[narma + (1L:ncxreg)], best_mod)
        else arimaSS(x, best_mod)

        best_sigma2 <- best_val[[1L]][1L]/n.used
        best_value <- 2 * n.used * best_res$value + n.used + n.used * log(2 * pi)
        all_values <- c(all_values, best_value)
      } else {  # New starting value

        # Many things can go wrong when using a random starting point as
        # initialization for ML estimation using stats::optim. Because the
        # CSS-ML approach of stats::arima is considered a reasonable approach,
        # we use that as a baseline and only get a different result if there is
        # no issues when performing the random restart algorithm. We check
        # ensure this using the try-catch statement below. We also suppress
        # warnings because each random starting point may produce a convergence
        # warning.

        suppressWarnings(
          restart_result <- tryCatch(
            {
              #### START

              if (include.mean) {
                new_init <- init
                new_init[mask] <- .sample_ARMA_coef(
                  arma = arma,
                  intercept = init[length(init)],
                  Mod_bounds = c(0.05, 0.95),
                  min_inv_root_dist = 0.01,
                  method = init_method
                )[mask]
              } else {
                new_init <- init
                new_init[mask] <- .sample_ARMA_coef(
                  arma = arma,
                  Mod_bounds = c(0.05, 0.95),
                  min_inv_root_dist = 0.01,
                  method = init_method
                )[mask]
              }

              if(transform.pars) {
                new_init <- .Call(C_ARIMA_Invtrans, new_init, arma)
                ## enforce invertibility
                if(arma[2L] > 0) {
                  ind <- arma[1L] + 1L:arma[2L]
                  new_init[ind] <- .maInvert(new_init[ind])
                }
                if(arma[4L] > 0) {
                  ind <- sum(arma[1L:3L]) + 1L:arma[4L]
                  new_init[ind] <- .maInvert(new_init[ind])
                }
              }
              trarma <- .Call(C_ARIMA_transPars, new_init, arma, transform.pars)
              mod <- stats::makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa, SSinit)
              res <- if(no.optim)
                list(convergence = 0, par = numeric(),
                     value = armafn(numeric(), as.logical(transform.pars)))
              else
                stats::optim(new_init[mask], armafn, method = optim.method,
                             hessian = TRUE, control = optim.control,
                             trans = as.logical(transform.pars))
              if(res$convergence > 0)
                warning(gettextf("possible convergence problem: optim gave code = %d",
                                 res$convergence), domain = NA)
              coef[mask] <- res$par
              if(transform.pars) {
                ## enforce invertibility
                if(arma[2L] > 0L) {
                  ind <- arma[1L] + 1L:arma[2L]
                  if(all(mask[ind]))
                    coef[ind] <- .maInvert(coef[ind])
                }
                if(arma[4L] > 0L) {
                  ind <- sum(arma[1L:3L]) + 1L:arma[4L]
                  if(all(mask[ind]))
                    coef[ind] <- .maInvert(coef[ind])
                }
                if(any(coef[mask] != res$par))  {  # need to re-fit
                  oldcode <- res$convergence
                  res <- stats::optim(coef[mask], armafn, method = optim.method,
                                      hessian = TRUE,
                                      control = list(maxit = 0L,
                                                     parscale = optim.control$parscale),
                                      trans = TRUE)
                  res$convergence <- oldcode
                  coef[mask] <- res$par
                }
                ## do it this way to ensure hessian was computed inside
                ## stationarity region
                A <- .Call(C_ARIMA_Gradtrans, as.double(coef), arma)
                A <- A[mask, mask]
                var <- crossprod(A, solve(res$hessian * n.used, A))
                coef <- .Call(C_ARIMA_undoPars, coef, arma)
              } else var <- if(no.optim) numeric() else solve(res$hessian * n.used)
              trarma <- .Call(C_ARIMA_transPars, coef, arma, FALSE)
              mod <- stats::makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa, SSinit)
              val <- if(ncxreg > 0L)
                arimaSS(x - xreg %*% coef[narma + (1L:ncxreg)], mod)
              else arimaSS(x, mod)
              sigma2 <- val[[1L]][1L]/n.used
              value <- 2 * n.used * res$value + n.used + n.used * log(2 * pi)

              list(
                i_trarma = trarma,
                i_mod = mod,
                i_val = val,
                i_sigma2 = sigma2,
                i_value = value,
                i_res = res,
                i_coef = coef,
                i_var = var,
                error = 0
              )
              #### End
            },
            error = function(e) list(i_trarma = trarma, i_mod = mod, i_val = NULL,
                                     i_sigma2 = NULL, i_value = Inf, i_res = NULL,
                                     i_coef = NULL, i_var = NULL, error = 1)
          )  # End trycatch
        )

        # Make sure the best fit also has a proper covariant matrix for coefficients.
        suppressWarnings(
          valid_test <- !is.null(restart_result$i_var) && !any(is.nan(sqrt(diag(restart_result$i_var))))
        )

        if (do_min_dist_test && !is.null(restart_result$i_coef)) {  # Check if there are (nearly) canceling roots;
          if (arma[1L] > 0 && arma[2L] > 0) {  # If there are both AR and MA coefs
            tmp_ar_pars <- restart_result$i_coef[1L:arma[1L]]
            tmp_ma_pars <- restart_result$i_coef[arma[1L] + 1L:arma[2L]]

            inv_ma_roots  <- 1 / polyroot(c(1, tmp_ma_pars))
            inv_ar_roots  <- 1 / polyroot(c(1, -tmp_ar_pars))
            inv_root_dist <- min(Mod(outer(inv_ar_roots, inv_ma_roots, FUN = '-')))
            valid_test <- valid_test && (inv_root_dist > min_inv_root_dist)
          }

          if (arma[3L] > 0 && arma[4L] > 0) {  # If there are both seasonal MA and seasonal AR coefs
            tmp_sar_pars <- restart_result$i_coef[sum(arma[1L:2L]) + 1L:arma[3L]]
            tmp_sma_pars <- restart_result$i_coef[sum(arma[1L:3L]) + 1L:arma[4L]]

            inv_sma_roots  <- 1 / polyroot(c(1, tmp_sma_pars))
            inv_sar_roots  <- 1 / polyroot(c(1, -tmp_sar_pars))
            inv_sroot_dist <- min(Mod(outer(inv_sar_roots, inv_sma_roots, FUN = '-')))
            valid_test <- valid_test && (inv_sroot_dist > min_inv_root_dist)
          }
        }

        if (do_maxroot_test && !is.null(restart_result$i_coef)) {  # Make sure all roots are not near boundary
          if (arma[1L] > 0) {  # If there are both AR coefs
            tmp_ar_pars <- restart_result$i_coef[1L:arma[1L]]
            inv_ar_roots  <- 1 / polyroot(c(1, -tmp_ar_pars))
            valid_test <- valid_test && (max(Mod(inv_ar_roots)) < max_inv_root)
          }

          if (arma[2L] > 0) {  # If there are MA coefs
            tmp_ma_pars <- restart_result$i_coef[arma[1L] + 1L:arma[2L]]
            inv_ma_roots  <- 1 / polyroot(c(1, tmp_ma_pars))
            valid_test <- valid_test && (max(Mod(inv_ma_roots)) < max_inv_root)
          }

          if (arma[3L] > 0) {  # If there are both seasonal MA and seasonal AR coefs
            tmp_sar_pars <- restart_result$i_coef[sum(arma[1L:2L]) + 1L:arma[3L]]
            inv_sar_roots  <- 1 / polyroot(c(1, -tmp_sar_pars))
            valid_test <- valid_test && (max(Mod(inv_sar_roots)) < max_inv_root)
          }

          if (arma[4L] > 0) {
            tmp_sma_pars <- restart_result$i_coef[sum(arma[1L:3L]) + 1L:arma[4L]]
            inv_sma_roots  <- 1 / polyroot(c(1, tmp_sma_pars))
            valid_test <- valid_test && (max(Mod(inv_sma_roots)) < max_inv_root)
          }
        }

        if (restart_result$error == 0 && restart_result$i_value + (eps_tol * 2) < best_value && valid_test) {
          best_coef <- restart_result$i_coef
          best_res <- restart_result$i_res
          best_var <- restart_result$i_var
          best_trarma <- restart_result$i_trarma
          best_mod <- restart_result$i_mod
          best_val <- restart_result$i_val
          best_sigma2 <- restart_result$i_sigma2
          best_value <- restart_result$i_value
          num_repeat <- 1
        } else {
          num_repeat <- num_repeat + 1
        }

        all_values <- c(all_values, best_value)

        converged <- num_repeat >= max_repeats
        # converged <- pnorm(2 * 0.3 * sqrt(i_start)) - pnorm(-2 * 0.3 * sqrt(i_start)) - (1 - get_rho_hat(all_values, eps_tol)/i_start)^(i_start + num_repeat) >= 0.99

        if (converged) {
          break
        }

      }
    }

    coef <- best_coef
    res <- best_res
    var <- best_var
    trarma <- best_trarma
    mod <- best_mod
    val <- best_val
    sigma2 <- best_sigma2
    value <- best_value

  }  # End of ML Estimation

  # value <- 2 * n.used * res$value + n.used + n.used * log(2 * pi)
  aic <- if(method != "CSS") value + 2*sum(mask) + 2 else NA
  nm <- NULL
  if (arma[1L] > 0L) nm <- c(nm, paste0("ar", 1L:arma[1L]))
  if (arma[2L] > 0L) nm <- c(nm, paste0("ma", 1L:arma[2L]))
  if (arma[3L] > 0L) nm <- c(nm, paste0("sar", 1L:arma[3L]))
  if (arma[4L] > 0L) nm <- c(nm, paste0("sma", 1L:arma[4L]))
  if (ncxreg > 0L) {
    nm <- c(nm, cn)
    if(!orig.xreg) {
      ind <- narma + 1L:ncxreg
      coef[ind] <- S$v %*% coef[ind]
      A <- diag(narma + ncxreg)
      A[ind, ind] <- S$v
      A <- A[mask, mask]
      var <- A %*% var %*% t(A)
    }
  }
  names(coef) <- nm
  if(!no.optim) dimnames(var) <- list(nm[mask], nm[mask])
  resid <- val[[2L]]
  stats::tsp(resid) <- xtsp
  class(resid) <- "ts"
  structure(list(coef = coef, sigma2 = sigma2, var.coef = var, mask = mask,
                 loglik = -0.5 * value, aic = aic, arma = arma,
                 residuals = resid, call = match.call(), series = series,
                 code = res$convergence, n.cond = ncond, nobs = n.used,
                 model = mod, x = x, num_starts = i_start,
                 all_values = -0.5 * all_values),
            class = c("Arima2", "Arima"))
}

