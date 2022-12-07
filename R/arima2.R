#' ARIMA Modelling of Time Series
#'
#' Fit an ARIMA model to a univariate time series. This function builds on
#' the ARIMA model fitting approach used in [stats::arima()] by fitting
#' model parameters via a random restart algorithm.
#'
#' @param nrestart Number of random restarts to use in fitting the model. If
#'    nrestart = 0L, then the function is equivalent to [stats::arima()].
#' @param method Fitting method. The default (unless there are missing values)
#'    is to use conditional-sum-of-squares to find starting values, then
#'    maximum likelihood.
#' @inheritParams stats::arima
#' @inherit stats::arima return
#' @export
#'
#' @examples
#' arima2(LakeHuron, order = c(2, 0, 1), nrestart = 10)
arima2 <- function(x, order = c(0L, 0L, 0L),
                   seasonal = list(order = c(0L, 0L, 0L), period = NA),
                   xreg = NULL, include.mean = TRUE, nrestart = 10,
                   transform.pars = TRUE,
                   fixed = NULL, init = NULL,
                   method = c("CSS-ML", "ML"),
                   n.cond,
                   SSinit = c("Gardner1980", "Rossignol2011"),
                   optim.method = "BFGS", optim.control = list(),
                   kappa = 1000000) {

  #  This function is based on the arima function of the stats package
  #  of R. Below the copright statement of the arima function is reproduced.
  #
  #  File src/library/stats/R/arima.R
  #  Part of the R package, http://www.R-project.org
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
  # Revised:

  C_TSconv <- utils::getFromNamespace("C_TSconv", "stats")
  C_getQ0 <- utils::getFromNamespace("C_getQ0", "stats")
  C_getQ0bis <- utils::getFromNamespace("C_getQ0bis", "stats")
  C_ARIMA_Like <- utils::getFromNamespace("C_ARIMA_Like", "stats")
  C_ARIMA_transPars <- utils::getFromNamespace("C_ARIMA_transPars", "stats")
  C_ARIMA_CSS <- utils::getFromNamespace("C_ARIMA_CSS", "stats")
  C_ARIMA_Invtrans <- utils::getFromNamespace("C_ARIMA_Invtrans", "stats")
  C_ARIMA_Gradtrans <- utils::getFromNamespace("C_ARIMA_Gradtrans", "stats")
  C_ARIMA_undoPars <- utils::getFromNamespace("C_ARIMA_undoPars", "stats")

  "%+%" <- function(a, b) .Call(C_TSconv, a, b)
  SSinit <- match.arg(SSinit)
  SS.G <- SSinit == "Gardner1980"
  upARIMA <- function(mod, phi, theta) {
    p <- length(phi)
    q <- length(theta)
    mod$phi <- phi
    mod$theta <- theta
    r <- max(p, q + 1L)
    if (p > 0)
      mod$T[1L:p, 1L] <- phi
    if (r > 1L)
      mod$Pn[1L:r, 1L:r] <- if (SS.G)
        .Call(C_getQ0, phi, theta)
    else .Call(C_getQ0bis, phi, theta, tol = 0)
    else mod$Pn[1L, 1L] <- if (p > 0)
      1/(1 - phi^2)
    else 1
    mod$a[] <- 0
    mod
  }
  arimaSS <- function(y, mod) {
    .Call(C_ARIMA_Like, y, mod, 0L, TRUE)
  }
  armafn <- function(p, trans) {
    par <- coef
    par[mask] <- p
    trarma <- .Call(C_ARIMA_transPars, par, arma, trans)
    if (is.null(Z <- tryCatch(upARIMA(mod, trarma[[1L]],
                                      trarma[[2L]]), error = function(e) NULL)))
      return(.Machine$double.xmax)
    if (ncxreg > 0)
      x <- x - xreg %*% par[narma + (1L:ncxreg)]
    res <- .Call(C_ARIMA_Like, x, Z, 0L, FALSE)
    s2 <- res[1L]/res[3L]
    0.5 * (log(s2) + res[2L]/res[3L])
  }
  armaCSS <- function(p) {
    par <- as.double(fixed)
    par[mask] <- p
    trarma <- .Call(C_ARIMA_transPars, par, arma, FALSE)
    if (ncxreg > 0)
      x <- x - xreg %*% par[narma + (1L:ncxreg)]
    res <- .Call(C_ARIMA_CSS, x, arma, trarma[[1L]], trarma[[2L]],
                 as.integer(ncond), FALSE)
    0.5 * log(res)
  }

  maInvert <- function(ma) {
    q <- length(ma)
    q0 <- max(which(c(1, ma) != 0)) - 1L
    if (!q0)
      return(ma)
    roots <- polyroot(c(1, ma[1L:q0]))
    ind <- Mod(roots) < 1
    if (all(!ind))
      return(ma)
    if (q0 == 1)
      return(c(1/ma[1L], rep.int(0, q - q0)))
    roots[ind] <- 1/roots[ind]
    x <- 1
    for (r in roots) x <- c(x, 0) - c(0, x)/r
    c(Re(x[-1L]), rep.int(0, q - q0))
  }

  getInits <- function(arma, init, mask) {
    # Create matrix that is just the init object repeated, will be modified as
    # needed.
    out_init <- matrix(
      rep(init, nrestart),
      nrow = nrestart,
      ncol = length(init),
      byrow = TRUE
    )

    if (arma[1L] > 0L) {  # arma[1L] contains number of AR coefficients.

      ar_ind <- 1:arma[1L]  # Get column index of AR coefficients

      if (any(mask[ar_ind])) {  # If not all are fixed, then estimate them.
        # n_ar <- sum(mask[ar_ind])  # Number of AR terms to estimate
        which_est_cols <- ar_ind[mask[ar_ind]]  # Of the AR columns, which are estimated.

        # Randomly generate AR coefficients, save them to appropriate columns
        out_init[, which_est_cols] <- matrix(
          stats::runif(nrestart * length(which_est_cols), min = -1, max = 1),
          nrow = nrestart, ncol = length(which_est_cols)
        )

        # Check if randomly generated AR coefficients make causal model
        bad_ar <- which(!apply(out_init[, ar_ind, drop = FALSE], 1, .arCheck))

        # Keep generating parameters until all AR coefficients correspond to
        # a causal model
        while (length(bad_ar) > 0) {

          # Generate new AR coefficients where needed
          out_init[bad_ar, which_est_cols] <- matrix(
            stats::runif(length(bad_ar) * length(which_est_cols), min = -1, max = 1),
            nrow = length(bad_ar), ncol = length(which_est_cols)
          )

          # Check if AR coefficients are causal
          bad_ar <- which(!apply(out_init[, ar_ind, drop = FALSE], 1, .arCheck))
        }
      }
    }

    # ARMA[2L] corresponds to the number of MA coefficients.
    if (arma[2L] > 0L) {

      # Get the indices for MA coefficients
      ma_ind <- arma[1L] + 1L:arma[2L]

      if (any(mask[ma_ind])) {  # If not all are fixed, then estimate them.
        which_est_cols <- ma_ind[mask[ma_ind]]  # Of the MA columns, which are estimated.

        # Generate random starting points for MA coefficients
        out_init[, which_est_cols] <- matrix(
          stats::runif(nrestart * length(which_est_cols), min = -1, max = 1),
          nrow = nrestart, ncol = length(which_est_cols)
        )

        out_init[, ma_ind] <- t(apply(out_init[, ma_ind, drop = FALSE], 1, maInvert))
      }
    }

    if (arma[3L] > 0L) {
      ar_seas_ind <- sum(arma[1L:2L]) + 1L:arma[3L]

      if (any(mask[ar_seas_ind])) {  # If not all are fixed, then estimate them.
        which_est_cols <- ar_seas_ind[mask[ar_seas_ind]]  # Of the AR columns, which are estimated.

        # Randomly generate AR coefficients, save them to appropriate columns
        out_init[, which_est_cols] <- matrix(
          stats::runif(nrestart * length(which_est_cols), min = -1, max = 1),
          nrow = nrestart, ncol = length(which_est_cols)
        )

        # Check if randomly generated AR coefficients correspond to causal model
        bad_ar <- which(!apply(out_init[, ar_seas_ind, drop = FALSE], 1, .arCheck))

        # Keep generating parameters until all AR coefficients correspond to
        # a causal model
        while (length(bad_ar) > 0) {

          # Generate new AR coefficients where needed
          out_init[bad_ar, which_est_cols] <- matrix(
            stats::runif(length(bad_ar) * length(which_est_cols), min = -1, max = 1),
            nrow = length(bad_ar), ncol = length(which_est_cols)
          )

          # Check if AR coefficients are causal
          bad_ar <- which(!apply(out_init[, ar_seas_ind, drop = FALSE], 1, .arCheck))
        }
      }
    }

    if (arma[4L] > 0L) {
      ma_seas_ind <- sum(arma[1L:3L]) + 1L:arma[4L]

      if (any(mask[ma_seas_ind])) {  # If not all are fixed, then randomly initialize them.
        which_est_cols <- ma_seas_ind[mask[ma_seas_ind]]  # Of the MA columns, which are estimated.

        # Generate random starting points for MA coefficients
        out_init[, which_est_cols] <- matrix(
          stats::runif(nrestart * length(which_est_cols), min = -1, max = 1),
          nrow = nrestart, ncol = length(which_est_cols)
        )

        out_init[, ma_seas_ind] <- t(apply(out_init[, ma_seas_ind, drop = FALSE], 1, maInvert))
      }
    }

    # Check if the model has an intercept, and if the intercept is estimated.
    if (include.mean && mask[length(mask)]) {

      # If intercept is fit, stay close to the intercept of init, since it's
      # unlikely to change much.
      out_init[, ncol(out_init)] <- out_init[, ncol(out_init)] + stats::rnorm(nrestart, sd = 0.05)
    }

    rbind(init, out_init)  # return the random initial values.
  }

  series <- deparse1(substitute(x))

  if (NCOL(x) > 1L) {
    stop("only implemented for univariate time series")
  }

  if (!is.null(fixed)) {
    warning("The random restart algorithm is not yet implemented for ARIMA
            models with fixed components. Setting nrester to 0L")
    nrestart <- 0L
  }

  method <- match.arg(method)
  x <- stats::as.ts(x)
  if (!is.numeric(x)) {
    stop("'x' must be numeric")
  }

  storage.mode(x) <- "double"
  dim(x) <- NULL
  n <- length(x)

  if (!missing(order)) {
    if (!is.numeric(order) || length(order) != 3L || any(order < 0)) {
      stop("'order' must be a non-negative numeric vector of length 3")
    }
  }
  if (!missing(seasonal))
    if (is.list(seasonal)) {
      if (is.null(seasonal$order))
        stop("'seasonal' must be a list with component 'order'")
      if (!is.numeric(seasonal$order) || length(seasonal$order) !=
          3L || any(seasonal$order < 0L))
        stop("'seasonal$order' must be a non-negative numeric vector of length 3")
    }
  else if (is.numeric(order)) {
    if (length(order) == 3L)
      seasonal <- list(order = seasonal)
    else ("'seasonal' is of the wrong length")
  }
  else stop("'seasonal' must be a list with component 'order'")
  if (is.null(seasonal$period) || is.na(seasonal$period) ||
      seasonal$period == 0)
    seasonal$period <- stats::frequency(x)
  arma <- as.integer(c(order[-2L], seasonal$order[-2L], seasonal$period,
                       order[2L], seasonal$order[2L]))
  narma <- sum(arma[1L:4L])
  xtsp <- stats::tsp(x)
  stats::tsp(x) <- NULL
  Delta <- 1
  for (i in seq_len(order[2L])) Delta <- Delta %+% c(1, -1)
  for (i in seq_len(seasonal$order[2L])) Delta <- Delta %+%
    c(1, rep.int(0, seasonal$period - 1), -1)
  Delta <- -Delta[-1L]
  nd <- order[2L] + seasonal$order[2L]
  n.used <- sum(!is.na(x)) - length(Delta)
  if (is.null(xreg)) {
    ncxreg <- 0L
  } else {
    nmxreg <- deparse1(substitute(xreg))
    if (NROW(xreg) != n)
      stop("lengths of 'x' and 'xreg' do not match")
    ncxreg <- NCOL(xreg)
    xreg <- as.matrix(xreg)
    storage.mode(xreg) <- "double"
  }
  class(xreg) <- NULL
  if (ncxreg > 0L && is.null(colnames(xreg)))
    colnames(xreg) <- if (ncxreg == 1L)
      nmxreg
  else paste0(nmxreg, 1L:ncxreg)
  if (include.mean && (nd == 0L)) {
    xreg <- cbind(intercept = rep(1, n), xreg = xreg)
    ncxreg <- ncxreg + 1L
  }
  if (method == "CSS-ML") {
    anyna <- anyNA(x)
    if (ncxreg)
      anyna <- anyna || anyNA(xreg)
    if (anyna)
      method <- "ML"

    ncond <- order[2L] + seasonal$order[2L] * seasonal$period
    ncond1 <- order[1L] + seasonal$period * seasonal$order[1L]
    ncond <- ncond + if (!missing(n.cond))
      max(n.cond, ncond1)
    else ncond1

  } else ncond <- 0

  if (is.null(fixed))
    fixed <- rep(NA_real_, narma + ncxreg)
  else if (length(fixed) != narma + ncxreg)
    stop("wrong length for 'fixed'")
  mask <- is.na(fixed)
  no.optim <- !any(mask)
  if (no.optim)
    transform.pars <- FALSE
  if (transform.pars) {
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
    if (order[2L] > 0L) {
      dx <- diff(dx, 1L, order[2L])
      dxreg <- diff(dxreg, 1L, order[2L])
    }
    if (seasonal$period > 1L && seasonal$order[2L] > 0) {
      dx <- diff(dx, seasonal$period, seasonal$order[2L])
      dxreg <- diff(dxreg, seasonal$period, seasonal$order[2L])
    }
    fit <- if (length(dx) > ncol(dxreg))
      stats::lm(dx ~ dxreg - 1, na.action = stats::na.omit)
    else list(rank = 0L)
    if (fit$rank == 0L) {
      fit <- stats::lm(x ~ xreg - 1, na.action = stats::na.omit)
    }
    isna <- is.na(x) | apply(xreg, 1L, anyNA)
    n.used <- sum(!isna) - length(Delta)
    init0 <- c(init0, coef(fit))
    ses <- summary(fit)$coefficients[, 2L]
    parscale <- c(parscale, 10 * ses)
  }
  if (n.used <= 0)
    stop("too few non-missing observations")
  if (!is.null(init)) {
    if (length(init) != length(init0))
      stop("'init' is of the wrong length")
    if (any(ind <- is.na(init)))
      init[ind] <- init0[ind]
    if (method == "ML") {
      if (arma[1L] > 0)
        if (!.arCheck(init[1L:arma[1L]]))
          stop("non-stationary AR part")
      if (arma[3L] > 0)
        if (!.arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
          stop("non-stationary seasonal AR part")
      if (transform.pars)
        init <- .Call(C_ARIMA_Invtrans, as.double(init),
                      arma)
    }
  }
  else init <- init0
  coef <- as.double(fixed)
  if (!("parscale" %in% names(optim.control)))
    optim.control$parscale <- parscale[mask]

  #### START

  if (method == "CSS-ML") {
    res <- if (no.optim)
      list(convergence = 0L, par = numeric(), value = armaCSS(numeric()))
    else {
      stats::optim(
        init[mask], armaCSS, method = optim.method,
        hessian = FALSE, control = optim.control
      )
    }

    if (res$convergence == 0)
      init[mask] <- res$par
    if (arma[1L] > 0)
      if (!.arCheck(init[1L:arma[1L]]))
        stop("non-stationary AR part from CSS")
    if (arma[3L] > 0)
      if (!.arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
        stop("non-stationary seasonal AR part from CSS")
    ncond <- 0L
  }
  if (transform.pars) {
    init <- .Call(C_ARIMA_Invtrans, init, arma)
    if (arma[2L] > 0) {
      ind <- arma[1L] + 1L:arma[2L]
      init[ind] <- maInvert(init[ind])
    }
    if (arma[4L] > 0) {
      ind <- sum(arma[1L:3L]) + 1L:arma[4L]
      init[ind] <- maInvert(init[ind])
    }
  }
  trarma <- .Call(C_ARIMA_transPars, init, arma, transform.pars)
  mod <- stats::makeARIMA(trarma[[1L]], trarma[[2L]], Delta,
                          kappa, SSinit)

  if (no.optim) {
    res <- list(
      convergence = 0, par = numeric(),
      value = armafn(numeric(),as.logical(transform.pars))
    )
  } else {

    if (nrestart != 0L) {

      best_val <- Inf
      coef_orig <- coef
      inits <- matrix(nrow = nrestart, ncol = length(init))
      Errs <- c()
      Err_MSG <- list()

      # Random Restart Algorithm.
      for (i in 1:nrestart) {

        if (i == 1L) {
          new_init <- init
        } else if (include.mean) {
          new_init <- .sample_ARMA_coef(arma, init[length(init)])
        } else {
          new_init <- .sample_ARMA_coef(arma)
        }

        inits[i, ] <- new_init

        suppressWarnings(
          res_temp <- tryCatch(
            list(fit = stats::optim(
              new_init, armafn, method = optim.method,
              hessian = TRUE, control = optim.control,
              trans = as.logical(transform.pars)
            ), e = 0),
            error = function(e) list(fit = list(value = Inf), e = 1, mes = e)
            # warning = function(w) list(fit = list(value = Inf), w = 1, e = 0)
          )
        )

        if (res_temp$e == 1) {
          Errs <- c(Errs, i)
          Err_MSG <- c(Err_MSG, res_temp$mes)
        }

        # If there were no issues with fitting from the starting point,
        # then check if it results in a better fit.
        if (res_temp$e != 1 && length(res_temp$fit$par) > 0) {
          coef_temp <- coef_orig
          coef_temp[mask] <- res_temp$fit$par
          if (transform.pars) {
            if (arma[2L] > 0L) {
              ind <- arma[1L] + 1L:arma[2L]
              if (all(mask[ind]))
                coef_temp[ind] <- maInvert(coef_temp[ind])
            }
            if (arma[4L] > 0L) {
              ind <- sum(arma[1L:3L]) + 1L:arma[4L]
              if (all(mask[ind]))
                coef_temp[ind] <- maInvert(coef_temp[ind])
            }
            if (any(coef_temp[mask] != res_temp$fit$par)) {
              oldcode <- res_temp$fit$convergence
              res_temp <- tryCatch(
                list(fit = stats::optim(
                  coef_temp[mask], armafn, method = optim.method,
                  hessian = TRUE,
                  control = list(maxit = 0L, parscale = optim.control$parscale),
                  trans = TRUE
                ), e = 0),
                error = function(e) list(fit = list(value = Inf, par = numeric(length(coef_temp[mask]))), e = 1, mes = e)
              )

              res_temp$fit$convergence <- oldcode
              coef_temp[mask] <- res_temp$fit$par
            }
            A_temp <- .Call(C_ARIMA_Gradtrans, as.double(coef_temp), arma)
            A_temp <- A_temp[mask, mask]
            var_temp <- tryCatch(
              crossprod(A_temp, solve(res_temp$fit$hessian * n.used, A_temp)),
              error = function(e) numeric()
            )
            coef_temp <- .Call(C_ARIMA_undoPars, coef_temp, arma)
          }
          else var_temp <- if (no.optim)
            numeric()
          else tryCatch(solve(res_temp$fit$hessian * n.used), error = function(e) numeric())
          trarma_temp <- .Call(C_ARIMA_transPars, coef_temp, arma, FALSE)
          mod_temp <- stats::makeARIMA(trarma_temp[[1L]], trarma_temp[[2L]], Delta,
                                       kappa, SSinit)
          val_temp <- if (ncxreg > 0L)
            arimaSS(x - xreg %*% coef_temp[narma + (1L:ncxreg)], mod_temp)
          else arimaSS(x, mod_temp)
          sigma2_temp <- val_temp[[1L]][1L]/n.used

          value_temp <- 2 * n.used * res_temp$fit$value + n.used + n.used * log(2 * pi)

          if (res_temp$e == 1 || length(var_temp) < 1) {
            Errs <- c(Errs, i)
            Err_MSG <- c(Err_MSG, res_temp$mes)
          }

          if (value_temp < best_val && length(var_temp) >= 1 && res_temp$e == 0) {
            best_val <- value_temp

            coef <- coef_temp
            res <- res_temp$fit
            A <- A_temp
            var <- var_temp
            trarma <- trarma_temp
            mod <- mod_temp
            val <- val_temp
            sigma2 <- sigma2_temp
            value <- value_temp
          }
        }
      }
    } else {  # nrestart == 0L

      # This is equivalent to stats::arima
      # Fit ARMA model with default initial values
      res <- stats::optim(
        init[mask], armafn, method = optim.method,
        hessian = TRUE, control = optim.control,
        trans = as.logical(transform.pars)
      )

      if (res$convergence > 0)
        warning(gettextf("possible convergence problem: optim gave code = %d",
                         res$convergence), domain = NA)
      coef[mask] <- res$par
      if (transform.pars) {
        if (arma[2L] > 0L) {
          ind <- arma[1L] + 1L:arma[2L]
          if (all(mask[ind]))
            coef[ind] <- maInvert(coef[ind])
        }
        if (arma[4L] > 0L) {
          ind <- sum(arma[1L:3L]) + 1L:arma[4L]
          if (all(mask[ind]))
            coef[ind] <- maInvert(coef[ind])
        }
        if (any(coef[mask] != res$par)) {
          oldcode <- res$convergence
          res <- stats::optim(
            coef[mask], armafn, method = optim.method,
            hessian = TRUE,
            control = list(maxit = 0L, parscale = optim.control$parscale),
            trans = TRUE
          )
          res$convergence <- oldcode
          coef[mask] <- res$par
        }
        A <- .Call(C_ARIMA_Gradtrans, as.double(coef), arma)
        A <- A[mask, mask]
        var <- crossprod(A, solve(res$hessian * n.used,
                                  A))
        coef <- .Call(C_ARIMA_undoPars, coef, arma)
      }
      else var <- if (no.optim)
        numeric()
      else solve(res$hessian * n.used)
      trarma <- .Call(C_ARIMA_transPars, coef, arma, FALSE)
      mod <- stats::makeARIMA(trarma[[1L]], trarma[[2L]], Delta,
                              kappa, SSinit)
      val <- if (ncxreg > 0L)
        arimaSS(x - xreg %*% coef[narma + (1L:ncxreg)], mod)
      else arimaSS(x, mod)
      sigma2 <- val[[1L]][1L]/n.used
    }
  }

  #### -end

  value <- 2 * n.used * res$value + n.used + n.used * log(2 * pi)
  aic <- value + 2 * sum(mask) + 2
  nm <- NULL
  if (arma[1L] > 0L)
    nm <- c(nm, paste0("ar", 1L:arma[1L]))
  if (arma[2L] > 0L)
    nm <- c(nm, paste0("ma", 1L:arma[2L]))
  if (arma[3L] > 0L)
    nm <- c(nm, paste0("sar", 1L:arma[3L]))
  if (arma[4L] > 0L)
    nm <- c(nm, paste0("sma", 1L:arma[4L]))
  if (ncxreg > 0L) {
    nm <- c(nm, cn)
    if (!orig.xreg) {
      ind <- narma + 1L:ncxreg
      coef[ind] <- S$v %*% coef[ind]
      A <- diag(narma + ncxreg)
      A[ind, ind] <- S$v
      A <- A[mask, mask]
      var <- A %*% var %*% t(A)
    }
  }
  names(coef) <- nm
  if (!no.optim) {
    dimnames(var) <- list(nm[mask], nm[mask])
  }
  resid <- val[[2L]]
  stats::tsp(resid) <- xtsp
  class(resid) <- "ts"
  structure(list(coef = coef, sigma2 = sigma2, var.coef = var,
                 mask = mask, loglik = -0.5 * value, aic = aic, arma = arma,
                 residuals = resid, call = match.call(), series = series,
                 code = res$convergence, n.cond = ncond, nobs = n.used,
                 model = mod, x = x, Errs = Errs, Msg = Err_MSG,
                 Inits = inits), class = c("Arima2", "Arima"))
}
