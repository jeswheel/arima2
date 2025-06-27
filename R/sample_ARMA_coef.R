#' Sample ARMA coef
#'
#' This function randomly samples the ARMA coefficients of a specified
#' ARMA model. The coefficients are sampled so that the resulting ARMA model
#' is both causal and invertible.
#'
#' For an ARMA model to be causal and invertible, the roots of the AR and MA
#' polynomials must lie outside the the complex unit circle. The AR and MA
#' polynomials are defined as:
#' \deqn{\phi(z) = 1 - \phi_1 z - \phi_2 z^2 - \ldots - \phi_p z^p}
#' \deqn{\theta(z) = 1 + \theta_1 z + \theta_2 z^2 + \ldots + \theta_q z^q}
#' where \eqn{z} is a complex number, \eqn{\phi_1, \ldots, \phi_p} are the
#' \eqn{p} AR coefficients of the ARMA model, and
#' \eqn{\theta_1, \ldots, \theta_p} are the \eqn{q} MA coefficients of the ARMA
#' model.
#'
#' This function implements two distinct sampling schemes.
#' \code{init_method = "DL"} will sample parameters using the Durbin-Levinson
#' algorithm, described by Monahan (1984). If \code{init_method = "UnifRoots"},
#' then inverted roots of AR and MA polynomials will be sampled uniformly from
#' the complex unit circle. In the later method, to ensure that the resulting
#' polynomial coefficients are real, we only sample half of the needed number
#' of complex roots, and set the remaining half to be the complex conjugate of
#' the sampled points. In the case where the number of coefficients is odd, the
#' remaining root is sampled uniformly, satisfying the `Mod_bounds` parameter.
#'
#' @param order A specification of the non-seasonal part of the ARIMA model:
#'    this is different than the `order` input of [stats::arima()], because
#'    the degree of differencing is not included. Thus order is a vector of
#'    length two of the form \eqn{(p, q)}.
#' @param seasonal A specification of the seasonal part of the ARIMA model. Can
#'    be either a vector of length 2, or a list with an `order` component; if
#'    a list, a period can be included, but it does not affect the function
#'    output.
#' @param n An integer indicating how many sets of ARMA coefficients should be
#'    sampled.
#' @param Mod_bounds Bounds on the magnitude of the roots.
#' @param min_inv_root_dist This parameter is included so as to help avoid ARMA
#'    models that contain parameter redundancy, if desired. Specifically,
#'    this parameter ensures that the minimum distance between any of the
#'    inverted roots in the AR and MA polynomials is greater than
#'    `min_inv_root_dist`. Inverted roots that are near each other leads to
#'    canceling or nearly canceling roots, effectively reducing the size of the
#'    ARMA model.
#' @param method  Method used to randomly sample parameter initializations.
#'    \code{init_method = "DL"} will sample parameters using the Durbin-Levinson
#'    algorithm, described by Monahan (1984). If
#'    \code{init_method = "UnifRoots"}, then inverted roots of AR and MA
#'    polynomials will be sampled uniformly from the complex unit circle.
#'
#' @returns a vector of randomly sampled ARMA coefficients.
#'
#' @references Monahan, John F. (1984) A note on enforcing stationarity in autoregressive-moving average models. \emph{Biometrika}, \bold{71}(2), 403--404.
#'
#' @export
#' @examples {
#' sample_ARMA_coef(
#'    order = c(2, 1),
#'    seasonal = list(order = c(1, 0), period = 2),
#'    n = 100
#' )
#' }
sample_ARMA_coef <- function(
    order = c(0L, 0L),
    seasonal = list(order = c(0L, 0L), period = NA),
    n = 1,
    Mod_bounds = c(0, 1),
    min_inv_root_dist = 0.0,
    method = c("UnifRoots", "DL")
) {

  if (!inherits(Mod_bounds, "numeric")) {
    stop("Argument 'Mod_bounds' must be a numeric vector of length 2.")
  } else if (length(Mod_bounds) != 2) {
    stop("Argument 'Mod_bounds' must be a numeric vector of length 2.")
  } else if (Mod_bounds[1] >= Mod_bounds[2]) {
    stop("First component of Mod_bounds must be less than the second component.")
  } else if (Mod_bounds[1] < 0) {
    stop("First component of Mod_bounds must be greater than or equal to zero.")
  } else if (Mod_bounds[2] > 1) {
    stop("Second component of Mod_bounds must be less than or equal to one.")
  }

  method <- match.arg(method)

  if (method == 'DL' && !all(Mod_bounds == c(0, 1))) {
    Mod_bounds <- c(0, 1)
    warning("If sampling method DL is used, Mod_bounds are set to c(0, 1).")
  }

  if (inherits(seasonal, "list")) {
    if (is.null(seasonal$order)) {
      stop("`seasonal` is missing component `order`.")
    } else if (length(seasonal$order) != 2) {
      stop("`seasonal$order` needs to be of length == 2")
    }

    seasonal <- seasonal$order
  }

  if (!is.numeric(min_inv_root_dist)) {
    stop("min_inv_root_dist must be numeric.")
  } else if (min_inv_root_dist < 0) {
    stop("min_inv_root_dist must be positive.")
  } else if (min_inv_root_dist > 1) {
    warning("min_inv_root_dist is too large. Consider lowering for improved efficiency.")
  }


  if (length(order) != 2) {
    stop("`order` should have length == 2.")
  }

  if (length(seasonal) != 2) {
    stop("Seasonal `order` should have length == 2.")
  }

  arma <- as.integer(
    c(order, seasonal, 1L, 0L, 0L)
  )

  if (n == 1L) {
    res <- .sample_ARMA_coef(arma, Mod_bounds = Mod_bounds, min_inv_root_dist = min_inv_root_dist, method = method)
  } else {
    res <- t(replicate(n, .sample_ARMA_coef(arma, Mod_bounds = Mod_bounds, min_inv_root_dist = min_inv_root_dist, method = method)))
  }

  res
}

#' Internal Sample ARMA Coef.
#'
#' @param arma vector of integers c(ar, ma, ar_seas, ma_seas, period, i, i_seas)
#' @param intercept If missing, the intercept is assumed fixed. Otherwise,
#'    new intercept values are sampled near the value of `intercept`.
#' @param min_inv_root_dist minimum allowable distance between inverted AR and
#'    MA polynomial roots, used to avoid models with parametric redundancies.
#' @param Mod_bounds Bounds on the magnitude of the roots.
#' @param method Type of sampling to be used. This mostly comes into play when
#'   doing the min_inv_root_dist test. If method is UnifRoot, then roots are
#'   sampled, checked, and then converted into coefficients. If method is
#'   DL, then parameters are simulated, then the roots are
#'   calculated and checked.
#'
#' This function does an internal check to make sure the resulting
#' coefficients do not result in AR and MA roots that approximately
#' cancel out. This is done by checking that the closest AR and MA
#' roots are at least 0.05 units apart.
#'
#' @return a vector of randomly sampled ARMA coefficients.
#' @noRd
.sample_ARMA_coef <- function(arma, intercept,
                              min_inv_root_dist,
                              Mod_bounds, method) {

  # get number of coefficients
  ar <- arma[1L]
  ma <- arma[2L]
  ar_seas <- arma[3L]
  ma_seas <- arma[4L]

  ar_coef <- numeric(ar)
  ma_coef <- numeric(ma)
  ar_seas_coef <- numeric(ar_seas)
  ma_seas_coef <- numeric(ma_seas)

  par_names <- c()

  if (ar != 0 & ma != 0) {  # Both AR and MA need sampling, check for canceling root.
    min_dist <- 0
    while (min_dist <= min_inv_root_dist) {

      if (method == 'DL') {
        # Using this method, we first sample coefs, then roots and calculate distance
        ar_coef <- .DLsample(n = ar, type = 'ar')
        names(ar_coef) <- paste0("ar", 1:ar)
        ma_coef <- .DLsample(n = ma, type = 'ma')
        names(ma_coef) <- paste0("ma", 1:ma)
        ar_inv_roots <- 1 / ARMApolyroots(ar_coef, type = 'AR')
        ma_inv_roots <- 1 / ARMApolyroots(ma_coef, type = 'MA')
        min_dist <- min(Mod(outer(ar_inv_roots, ma_inv_roots, FUN = '-')))
      } else {
        # Using other methods, we sample roots, calculate distance, then convert to coefs.
        ar_inv_roots <- .sample_inv_roots(n = ar, Mod_bounds = Mod_bounds)
        ma_inv_roots <- .sample_inv_roots(n = ma, Mod_bounds = Mod_bounds)
        min_dist <- min(Mod(outer(ar_inv_roots, ma_inv_roots, FUN = '-')))
      }
    }

    if (method != 'DL') {  # Need to get coefs if not "DL" method.
      ar_coef <- .roots2coef(ar_inv_roots, type = 'ar')
      ma_coef <- .roots2coef(ma_inv_roots, type = 'ma')
    }

    par_names <- c(par_names, paste0("ar", 1:arma[1L]))
    par_names <- c(par_names, paste0("ma", 1:arma[2L]))
  } else if (ar != 0) {
    par_names <- c(par_names, paste0("ar", 1:arma[1L]))

    if (method == "DL") {
      ar_coef <- .DLsample(ar, type = 'ar')
    } else {
      ar_coef <- .roots2coef(.sample_inv_roots(ar, Mod_bounds = Mod_bounds), type = 'ar')
    }
  } else if (ma != 0) {
    par_names <- c(par_names, paste0("ma", 1:arma[2L]))

    if (method == "DL") {
      ma_coef <- .DLsample(ma, type = 'ma')
    } else {
      ma_coef <- .roots2coef(.sample_inv_roots(ma, Mod_bounds = Mod_bounds), type = 'ma')
    }
  }

  if (ar_seas != 0 & ma_seas != 0) {  # Both AR and MA need sampling, check for canceling root.
    min_dist <- 0
    while (min_dist <= min_inv_root_dist) {

      if (method == 'DL') {
        ar_seas_coef <- .DLsample(n = ar_seas, type = 'ar')
        names(ar_seas_coef) <- paste0("ar", 1:ar_seas)
        ma_seas_coef <- .DLsample(n = ma_seas, type = 'ma')
        names(ma_seas_coef) <- paste0("ma", 1:ma_seas)
        ar_inv_roots <- 1 / ARMApolyroots(ar_seas_coef, type = 'AR')
        ma_inv_roots <- 1 / ARMApolyroots(ma_seas_coef, type = 'MA')
        min_dist <- min(Mod(outer(ar_inv_roots, ma_inv_roots, FUN = '-')))
      } else {
        # Resample everything. Condition is not often triggered, and no need to worry about ensuring complex conjugate pairs.
        ar_inv_roots <- .sample_inv_roots(n = ar_seas, Mod_bounds = Mod_bounds)
        ma_inv_roots <- .sample_inv_roots(n = ma_seas, Mod_bounds = Mod_bounds)
        min_dist <- min(Mod(outer(ar_inv_roots, ma_inv_roots, FUN = '-')))
      }
    }

    if (method != 'DL') {  # Need to get coefs if not "DL" method.
      ar_seas_coef <- .roots2coef(ar_inv_roots, type = 'ar')
      ma_seas_coef <- .roots2coef(ma_inv_roots, type = 'ma')
    }

    par_names <- c(par_names, paste0("ar_seas", 1:arma[3L]))
    par_names <- c(par_names, paste0("ma_seas", 1:arma[4L]))
  } else if (ar_seas != 0) {

    par_names <- c(par_names, paste0("ar_seas", 1:arma[3L]))

    if (method == "DL") {
      ar_seas_coef <- .DLsample(ar_seas, type = 'ar')
    } else {
      ar_seas_coef <- .roots2coef(.sample_inv_roots(ar_seas, Mod_bounds = Mod_bounds), type = 'ar')
    }

  } else if (ma_seas != 0) {

    if (method == "DL") {
      ma_seas_coef <- .DLsample(ma_seas, type = 'ma')
    } else {
      ma_seas_coef <- .roots2coef(.sample_inv_roots(ma_seas, Mod_bounds = Mod_bounds), type = 'ma')
    }

    par_names <- c(par_names, paste0("ma_seas", 1:arma[4L]))
  }

  out <- if (!missing(intercept)) {
    c(ar_coef, ma_coef, ar_seas_coef, ma_seas_coef, stats::rnorm(1, intercept, 0.05))
  } else {
    c(ar_coef, ma_coef, ar_seas_coef, ma_seas_coef)
  }

  if (!missing(intercept)) {
    par_names <- c(par_names, "intercept")
  }

  names(out) <- par_names

  return(out)
}

#' Roots 2 coef
#'
#' This function converts the inverse roots of an AR or MA polynomial and
#' calculates the corresponding polynomial, i.e., the coefficients of the
#' ARMA model
#'
#' @param inv_roots vector of inverse polynomial roots. Each should be complex
#'    numbers that lie within the complex-unit circle.
#' @param type string of "ar" or "ma", indicates whether or not the polynomial
#'    coefficients are negative or positive.
#'
#' @return Coefficients of an AR or MA model that corresponds to the inverse
#'    roots.
#'
#' @examples arima2:::.roots2coef(.sample_inv_roots(2))
#' @noRd
.roots2coef <- function(inv_roots, type = "ar") {

  if (!type %in% c("ar", "ma")) {
    stop("Type must be 'ar' or 'ma'")
  }

  x <- 1
  roots <- 1 / inv_roots
  for (r in roots) x <- c(x, 0) - c(0, x)/r
  coefs <- Re(x[-1L])

  if (type == 'ar') {
    -coefs
  } else {
    coefs
  }
}

#' Sample Inverse Roots
#'
#' This function samples inverse roots of an AR or MA polynomial so that the
#' reciprocal of the root lies within the complex unit circle.
#'
#' The sampling is either done by sampling a radius from a uniform(0, 1)
#' or a beta(2.5, 4) distribution. We chose beta(2.5, 4), because inverse roots
#' with radius close to 1 result in potentially problematic ARMA coefficients,
#' and values close to zero result in ARMA coefficients that are very similar
#' and near zero, which defeats the purpose of the random restart algorithm.
#'
#' Edit: Jun 26, 2025. The new coefficient sampling algorithm (Durbin-Levinson)
#' does not require sampling roots. Instead
#'
#' @param n number of roots to sample (equal to the number of coefficients)
#' @param type string. If type == "beta", the radius is sampled with a beta
#'    distribution. Otherwise, uniform sampling is used.
#'
#' @return random inverse roots that lie inside the complex unit circle
#' @noRd
#'
#' @examples arima2:::.sample_inv_roots(2, type = 'unif', Mod_bounds = c(0.05, 0.95))
.sample_inv_roots <- function(n, type = 'unif', Mod_bounds = c(0, 1)) {

  n_complex_pairs <- stats::rbinom(1, size = n %/% 2, prob = 1 - sqrt(0.5)) # 1 - sqrt(0.5)
  n_real_pairs <- max(0, n %/% 2 - n_complex_pairs)

  if (type == 'beta') {
    R <- stats::rbeta(n_complex_pairs, 2.5, 4)
    Theta <- pi * stats::runif(n_complex_pairs)
  } else if (type == 'unif') {
    R <- stats::runif(n_complex_pairs, Mod_bounds[1], Mod_bounds[2])
    Theta <- pi * stats::runif(n_complex_pairs)
  } else {
    stop("Not valid sampling type")
  }

  complex_roots <- complex(real = R * cos(Theta), imaginary = R * sin(Theta))
  inv_roots <- c(complex_roots, Conj(complex_roots))

  real_pairs <- matrix(nrow = n_real_pairs, ncol = 2)

  same_sign <- stats::rbinom(n_real_pairs, 1, 2 * sqrt(0.5) - 1) |> as.logical()
  flip_col <- sample(c(1, 2), size = sum(same_sign), replace = TRUE)
  real_pairs[same_sign, ] <- stats::runif(2*sum(same_sign), Mod_bounds[1], Mod_bounds[2]) * sample(c(-1, 1), size = sum(same_sign), replace = TRUE)
  # real_pairs[same_sign, ] <- real_pairs[same_sign, ] * sample(c(-1, 1), size = sum(same_sign), replace = TRUE)
  real_pairs[!same_sign, ] <- sample(c(-1, 1), 2 * (n_real_pairs - sum(same_sign)), replace = TRUE) * stats::runif(2 * (n_real_pairs - sum(same_sign)), Mod_bounds[1], Mod_bounds[2])
  real_pairs[cbind(which(same_sign), flip_col)] <- -1 * real_pairs[cbind(which(same_sign), flip_col)]
  inv_roots <- c(inv_roots, as.numeric(real_pairs))

  if (n %% 2 == 1) {
    re_inv_root <- sample(c(-1, 1), 1) * stats::runif(1, min = Mod_bounds[1], Mod_bounds[2])
    inv_roots <- c(inv_roots, re_inv_root)
  }

  inv_roots
}

#' Sample AR / MA coefficients using Durbin-Levinson
#'
#' @param n Number of coefficients to sample
#' @param type Character vector, of type c('ar', 'ma')
#' @param gamma bounds on the sampled PACF
#' @param theta bounds on the sampled PACF
#'
#' @returns A vector of length n containing AR or MA coefficients, specified by
#'    type
#' @noRd
#'
#' @examples arima2:::.DLsample(3, type = 'ar')
.DLsample <- function(n, type, gamma = 0.01, theta = 0.005) {
  phi_kk <- numeric(n)
  phi_kk[1] <- stats::runif(1, min = -1 + gamma, max = 1 - gamma)
  if (n > 1) phi_kk[2:n] <- stats::runif(n-1, min = -1+gamma+theta, max = 1 - gamma - theta)

  tmp <- .Call(C_ARIMA_transPars, atanh(phi_kk), as.integer(c(n, 0, 0, 0, 0, 0, 0)), TRUE)

  if (type == 'ar') tmp[[1L]] else -tmp[[1L]]
}

#' AR Check
#'
#' A function to check whether the ARMA model is causal, by checking the AR
#' polynomial has roots that lie outside the unit circle.
#'
#' @param ar vector of coefficients of an ARMA model.
#' @param Mod_bounds Additional argument used for sampling coefficients.
#'    This ensures that the sampled coefficients result in a reasonably
#'    well-behaved ARMA model, where each coefficient has a non-negligible
#'    effect, and the model isn't approaching stationarity.
#'
#' @return Boolean. TRUE if AR coefficients correspond to a causal AR model.
#' @noRd
#'
#' @examples arima2:::.arCheck(c(0.94, -0.1), Mod_bounds = c(0, 1))
.arCheck <- function(ar, Mod_bounds = c(0, 1)) {
  p <- max(which(c(1, -ar) != 0)) - 1

  if (!p) {
    return(TRUE)
  }

  root_size <- Mod(polyroot(c(1, -ar[1L:p])))

  all(root_size > 1) && min(1 / root_size) > Mod_bounds[1] && max(1 / root_size) < Mod_bounds[2]
}

#' MA Check
#'
#' A function to check whether the ARMA model is invertible, by checking the MA
#' polynomial has roots that lie outside the unit circle.
#'
#' @param ma vector of coefficients of an ARMA model.
#' @param Mod_bounds Additional argument used for sampling coefficients.
#'    This ensures that the sampled coefficients result in a reasonably
#'    well-behaved ARMA model, where each coefficient has a non-negligible
#'    effect, and the model isn't approaching stationarity.
#'
#' @return Boolean. TRUE if MA coefficients correspond to an invertible MA model.
#' @noRd
#'
#' @examples arima2:::.maCheck(c(0.94, -0.1), Mod_bounds = c(0, 1))
.maCheck <- function(ma, Mod_bounds = c(0, 1)) {
  p <- max(which(c(1, ma) != 0)) - 1

  if (!p) {
    return(TRUE)
  }


  root_size <- Mod(polyroot(c(1, ma[1L:p])))

  all(root_size > 1) && min(1 / root_size) > Mod_bounds[1] && max(1 / root_size) < Mod_bounds[2]
}

#' MA Invert
#'
#' A function to invert the MA coefficients. This was copy and pasted from
#' stats::arima
#'
#' @param ma vector of coefficients of an ARMA model.
#' @return ma coefficients.
#' @noRd
.maInvert <- function(ma) {

  #  This function is based on the arima function of the stats package
  #  of R. Below the copright statement of the arima function is reproduced.
  #
  #  File src/library/stats/R/arima.R
  #  Part of the R package, https://www.R-project.org
  #
  #  Copyright (C) 2002-2015 The R Core Team
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
  #  https://www.R-project.org/Licenses/

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
