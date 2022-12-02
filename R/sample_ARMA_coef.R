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
#' ARMA coefficients are sampled by sampling inverse roots to be inside the
#' complex unit circle, and then calculating the resulting polynomial.
#' To ensure that the resulting polynomial coefficients are real, we only sample
#' half of the needed number of complex roots, and set the remaining half to be
#' the complex conjugate of the sampled points. In the case where the number of
#' coefficients is odd, the first coefficient is sampled as a real number
#' between 0-1, the remaining roots are sampled as described, and then the
#' resulting polynomial is checked to ensure all roots lie outside the unit
#' circle; if not, we resample the coefficients.
#'
#' @param order A specification of the non-seasonal part of the ARIMA model:
#'    this is different than the `order` input of [stats::arima()], because
#'    the degree of differencing is not included. Thus order is a vector of
#'    length two of the form \eqn{(p, q)}.
#' @param seasonal A specification of the seasonal part of the ARIMA model. Can
#'    be either a vector of length 2, or a list with an `order` component; if
#'    a list, aperiod can be included, but it does not affect the function
#'    output.
#' @param n An integer indicating how many sets of ARMA coefficients should be
#'    sampled.
#'
#' @return a vector of randomly sampled ARMA coefficients.
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
    n = 1
) {

  if (inherits(seasonal, "list")) {
    if (is.null(seasonal$order)) {
      stop("`seasonal` is missing component `order`.")
    } else if (length(seasonal$order) != 2) {
      stop("`seasonal$order` needs to be of length == 2")
    }

    seasonal <- seasonal$order
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

  res_names <- 0

  if (order[1L] != 0) {
    res_names <- paste0("ar", 1:order[1L])
  }

  if (order[2L] != 0) {
    res_names <- c(res_names, paste0("ma", 1:order[2L]))
  }

  if (seasonal[1L] != 0) {
    res_names <- c(res_names, paste0("ar_seas", 1:seasonal[1L]))
  }

  if (seasonal[2L] != 0) {
    res_names <- c(res_names, paste0("ma_seas", 1:seasonal[2L]))
  }

  if (n == 1L) {
    res <- .sample_ARMA_coef(arma)
    names(res) <- res_names
  } else {
    res <- t(replicate(n, .sample_ARMA_coef(arma)))
    colnames(res) <- res_names
  }

  res
}

#' Internal Sample ARMA Coef.
#'
#' @param arma vector of integers c(ar, ma, ar_seas, ma_seas, period, i, i_seas)
#' @param intercept If missing, the intercept is assumed fixed. Otherwise,
#'    new intercept values are sampled near the value of `intercept`.
#'
#' @return a vector of randomly sampled ARMA coefficients.
#' @noRd
.sample_ARMA_coef <- function(arma, intercept) {

  # get number of coefficients
  ar <- arma[1L]
  ma <- arma[2L]
  ar_seas <- arma[3L]
  ma_seas <- arma[4L]

  if (ar %% 2L == 0L) {  # Even number

    if (ar == 0L) {
      ar_coef <- c()
    } else {
      ar_inv_roots <- .sample_inv_roots(ar / 2L)
      ar_inv_roots_conj <- Conj(ar_inv_roots)

      ar_coef <- Re(.roots2poly(c(ar_inv_roots, ar_inv_roots_conj), type = 'ar'))
    }

  } else {  # Odd number

    ar_coef <- stats::runif(1L, min = -1, max = 1)

    if (ar > 1L) {

      ar_coef <- rep(1, ar)

      while(!.arCheck(ar_coef)) {
        ar_inv_roots <- .sample_inv_roots((ar - 1L) / 2L)
        ar_inv_roots_conj <- Conj(ar_inv_roots)

        ar_coef <- Re(
          .roots2poly(c(stats::runif(1L, -1, 1), ar_inv_roots, ar_inv_roots_conj), type = 'ar')
        )
      }
    }
  }

  if (ma %% 2L == 0L) {  # Even number

    if (ma == 0L) {
      ma_coef <- c()
    } else {
      ma_inv_roots <- .sample_inv_roots(ma / 2L)
      ma_inv_roots_conj <- Conj(ma_inv_roots)

      ma_coef <- Re(
        .roots2poly(c(ma_inv_roots, ma_inv_roots_conj), type = 'ma')
      )
    }

  } else {  # Odd number

    ma_coef <- stats::runif(1L, -1, 1)

    if (ma > 1L) {

      ma_coef <- rep(1, ma)

      while(!.maCheck(ma_coef)) {
        ma_inv_roots <- .sample_inv_roots((ma - 1L) / 2L)
        ma_inv_roots_conj <- Conj(ma_inv_roots)

        ma_coef <- Re(
          .roots2poly(
            c(stats::runif(1L, -1, 1), ma_inv_roots, ma_inv_roots_conj),
            type = 'ma'
          )
        )
      }
    }
  }


  if (ar_seas %% 2L == 0L) {  # Even number

    if (ar_seas == 0L) {
      ar_seas_coef <- c()
    } else {
      ar_seas_inv_roots <- .sample_inv_roots(ar_seas / 2L)
      ar_seas_inv_roots_conj <- Conj(ar_seas_inv_roots)

      ar_seas_coef <- Re(
        .roots2poly(c(ar_seas_inv_roots, ar_seas_inv_roots_conj), type = 'ar')
      )
    }

  } else {  # Odd number
    ar_seas_coef <- stats::runif(1L, -1, 1)

    if (ar_seas > 1L) {

      ar_seas_coef <- rep(1, ar_seas)

      while(!.arCheck(ar_seas_coef)) {
        ar_seas_inv_roots <- .sample_inv_roots((ar_seas - 1L) / 2L)
        ar_seas_inv_roots_conj <- Conj(ar_seas_inv_roots)

        ar_seas_coef <- Re(
          .roots2poly(
            c(stats::runif(1L, -1, 1), ar_seas_inv_roots, ar_seas_inv_roots_conj),
            type = 'ar')
        )

        if(max(abs(Im(ar_seas_coef))) > 1e-12) {
          warning("Check AR coefficients.")
        }

        ar_seas_coef <- Re(ar_seas_coef)
      }
    }
  }

  if (ma_seas %% 2L == 0L) {  # Even number

    if (ma_seas == 0L) {
      ma_seas_coef <- c()
    } else {
      ma_seas_inv_roots <- .sample_inv_roots(ma_seas / 2L)
      ma_seas_inv_roots_conj <- Conj(ma_seas_inv_roots)

      ma_seas_coef <- Re(
        .roots2poly(c(ma_seas_inv_roots, ma_seas_inv_roots_conj), type = 'ma')
      )
    }

  } else {  # Odd number

    ma_seas_coef <- stats::runif(1L, -1, 1)

    if (ma_seas > 1L) {

      ma_seas_coef <- rep(1, ma_seas)

      while(!.maCheck(ma_seas_coef)) {
        ma_seas_inv_roots <- .sample_inv_roots((ma_seas - 1L) / 2L)
        ma_seas_inv_roots_conj <- Conj(ma_seas_inv_roots)

        ma_seas_coef <- Re(
          .roots2poly(
            c(stats::runif(1L, -1, 1), ma_seas_inv_roots, ma_seas_inv_roots_conj),
            type = 'ma'
          )
        )
      }
    }
  }

  if (!missing(intercept)) {
    c(ar_coef, ma_coef, ar_seas_coef, ma_seas_coef, stats::rnorm(1, intercept, 0.05))
  } else {
    c(ar_coef, ma_coef, ar_seas_coef, ma_seas_coef)
  }
}

#' Roots 2 poly
#'
#' This function converts the inverse roots of an AR or MA polynomial and
#' calculates the corresponding polynomial, i.e., the coefficients of the
#' ARMA model
#'
#'
#' @param inv_roots vector of inverse polynomial roots. Each should be complex
#'    numbers that lie within the complex-unit circle.
#' @param type string of "ar" or "ma", indicates whether or not the polynomial
#'    coefficients are negative or positive.
#'
#' @return Coefficients of an AR or MA model that corresponds to the inverse
#'    roots.
#'
#' @examples .roots2poly(.sample_inv_roots(2))
#' @noRd
.roots2poly <- function(inv_roots, type = "ar") {

  if (!type %in% c("ar", "ma")) {
    stop("Type must be 'ar' or 'ma'")
  }

  n_coef <- length(inv_roots)
  coefs <- numeric(n_coef)
  for (i in 1:n_coef) {
    coefs[i] <- utils::combn(x = inv_roots, m = i) |>
      apply(2, prod) |>
      sum() * (-1)^i
  }

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
#' @param n number of roots to sample (equal to the number of coefficients)
#' @param type string. If type == "beta", the radius is sampled with a beta
#'    distribution. Otherwise, uniform sampling is used.
#'
#' @return random inverse roots that lie inside the complex unit circle
#' @noRd
#'
#' @examples .sample_inv_roots(2)
.sample_inv_roots <- function(n, type = 'beta') {

  if (type == 'beta') {
    R <- stats::rbeta(n, 2.5, 4)
    Theta <- pi * stats::runif(n)
  } else {
    R <- stats::runif(n)
    Theta <- pi * stats::runif(n)
  }

  complex(real = R * cos(Theta), imaginary = R * sin(Theta))
}

#' AR Check
#'
#' A function to check whether the ARMA model is causal, by checking the AR
#' polynomial has roots that lie outside the unit circle.
#'
#' @param ar vector of coefficients of an ARMA model.
#'
#' @return Boolean. TRUE if AR coefficients correspond to a causal AR model.
#' @noRd
#'
#' @examples .arCheck(c(0.94, -0.1))
.arCheck <- function(ar) {
  p <- max(which(c(1, -ar) != 0)) - 1

  if (!p) {
    return(TRUE)
  }

  all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
}

#' MA Check
#'
#' A function to check whether the ARMA model is invertible, by checking the MA
#' polynomial has roots that lie outside the unit circle.
#'
#' @param ma vector of coefficients of an ARMA model.
#'
#' @return Boolean. TRUE if MA coefficients correspond to an invertible MA model.
#' @noRd
#'
#' @examples .maCheck(c(0.94, -0.1))
.maCheck <- function(ma) {
  p <- max(which(c(1, ma) != 0)) - 1

  if (!p) {
    return(TRUE)
  }

  all(Mod(polyroot(c(1, ma[1L:p]))) > 1)
}

