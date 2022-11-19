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
#'
#' @param arma vector of integers c(ar, ma, ar_seas, ma_seas, period, i, i_seas)
#' @param intercept If missing, the intercept is assumed fixed. Otherwise,
#'    new intercept values are sampled near the value of `intercept`.
#'
#' @return a vector of randomly sampled ARMA coefficients.
#'
#' @noRd
#' @examples sample_ARMA_coef(c(2, 1, 1, 0, 1, 0, 0), 153.1)
sample_ARMA_coef <- function(arma, intercept) {

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

    ar_coef <- runif(1L)

    if (ar > 1L) {

      ar_coef <- rep(1, ar)

      while(!.arCheck(ar_coef)) {
        ar_inv_roots <- .sample_inv_roots((ar - 1L) / 2L)
        ar_inv_roots_conj <- Conj(ar_inv_roots)

        ar_coef <- Re(
          .roots2poly(c(runif(1), ar_inv_roots, ar_inv_roots_conj), type = 'ar')
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

    ma_coef <- runif(1L)

    if (ma > 1L) {

      ma_coef <- rep(1, ma)

      while(!.maCheck(ma_coef)) {
        ma_inv_roots <- .sample_inv_roots((ma - 1L) / 2L)
        ma_inv_roots_conj <- Conj(ma_inv_roots)

        ma_coef <- Re(
          .roots2poly(
            c(runif(1), ma_inv_roots, ma_inv_roots_conj),
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
    ar_seas_coef <- runif(1L)

    if (ar_seas > 1L) {

      ar_seas_coef <- rep(1, ar_seas)

      while(!.arCheck(ar_seas_coef)) {
        ar_seas_inv_roots <- .sample_inv_roots((ar_seas - 1L) / 2L)
        ar_seas_inv_roots_conj <- Conj(ar_seas_inv_roots)

        ar_seas_coef <- Re(
          .roots2poly(
            c(runif(1), ar_seas_inv_roots, ar_seas_inv_roots_conj),
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

    ma_seas_coef <- runif(1L)

    if (ma_seas > 1L) {

      ma_seas_coef <- rep(1, ma_seas)

      while(!.maCheck(ma_seas_coef)) {
        ma_seas_inv_roots <- .sample_inv_roots((ma_seas - 1L) / 2L)
        ma_seas_inv_roots_conj <- Conj(ma_seas_inv_roots)

        ma_seas_coef <- Re(
          .roots2poly(
            c(runif(1), ma_seas_inv_roots, ma_seas_inv_roots_conj),
            type = 'ma'
          )
        )
      }
    }
  }

  if (!missing(intercept)) {
    c(ar_coef, ma_coef, ar_seas_coef, ma_seas_coef, rnorm(1, intercept, 0.05))
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
    R <- rbeta(n, 2.5, 4)
    Theta <- pi * runif(n)
  } else {
    R <- runif(n)
    Theta <- pi * runif(n)
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

