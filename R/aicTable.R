#' ARIMA AIC table
#'
#' Construct table of AIC for all combinations 0<=p<=P and 0<=q<=Q
#'
#' This function creates an AIC table for ARMA models of varying sizes.
#' Each row for the table corresponds to a different AR value, and each column
#' of the table corresponds to a different MA value.
#'
#' @param data a time series object, or a dataset that can be used as input into
#'    the [arima] function.
#' @param P a positive integer value representing the maximum number of AR
#'    coefficients that should be included in the table.
#' @param Q a positive integer value representing the maximum number of MA
#'    coefficients that should be included in the table.
#' @param D a positive integer value representing the degree of differencing
#' @param ic Information criterion to be used in the table.
#' @param ... Additional arguments passed to [arima()].
#'
#' @returns A matrix containing the model AIC values.
#' @export
#' @examples
#' set.seed(654321)
#' aicTable(presidents, 3, 2)
aicTable <- function(data, P, Q, D = 0, ic = c('aic', 'aicc'), ...){

  ic <- match.arg(ic)

  if (!is.numeric(P) | !is.numeric(Q) | !is.numeric(D)) {
    stop("'P', 'Q' and 'D' must be numeric.")
  }

  P <- as.integer(P)
  Q <- as.integer(Q)
  D <- as.integer(D)

  table <- matrix(NA, (P + 1), (Q + 1))
  for (p in 0:P) {
    for (q in 0:Q) {
      mod <- arima(data, order = c(p, D, q), ...)

      val <- mod$aic

      if (ic == "aicc") {
        k <- sum(mod$mask) + 1
        val <- val + (2 * k^2 + 2 * k) / (mod$nobs - k - 1)
      }

      table[p + 1, q + 1] <- val
    }
  }
  dimnames(table) <- list(paste("AR", 0:P, sep = ""), paste("MA", 0:Q, sep = ""))
  table
}

#' Check Table
#'
#' This function is used to check the consistency of an AIC table generated
#' using the aicTable function (above).
#'
#' This function was primarily implemented to help with the ArXiv paper that
#' describes this package, and for that reason aicc check isn't implemented.
#'
#' @param data a time series object, or a dataset that can be used as input into
#'    the [arima] function.
#' @param P a positive integer value representing the maximum number of AR
#'    coefficients that should be included in the table.
#' @param Q a positive integer value representing the maximum number of MA
#'    coefficients that should be included in the table.
#' @param D a positive integer value representing the degree of differencing
#' @param method string that must be "arima2" or "stats", indicating which
#'    package should be used to fit the ARIMA model.
#' @param eps_tol Tolerance for accepting a new solution to be better than a
#'    previous solution in terms of log-likelihood. The default corresponds to a
#'    one ten-thousandth unit increase in log-likelihood.
#' @param ... Additional arguments passed to either [stats::arima()] or
#'    [arima()], depending on which method is called.
#'
#' @return Boolean. True if the table is consistent in the sense that larger
#'    models have likelihood that is greater than or equal to all smaller
#'    models.
#' @noRd
#'
#' @examples arima2:::.checkTable(presidents, 3, 2)
.checkTable <- function(data, P, Q, D = 0, method = 'arima2', eps_tol = 1e-4, ...) {

  is_consistent = TRUE

  if (!is.numeric(P) | !is.numeric(Q) | !is.numeric(D)) {
    stop("'P', 'Q' and 'D' must be numeric.")
  }

  P <- as.integer(P)
  Q <- as.integer(Q)
  D <- as.integer(D)

  table <- matrix(NA, (P + 1), (Q + 1))
  for(p in 0:P) {
    for(q in 0:Q) {

      if (method == 'arima2') {
        table[p + 1, q + 1] <- arima(data, order = c(p, D, q), ...)$loglik
      } else {
        table[p + 1, q + 1] <- stats::arima(data, order = c(p, D, q), ...)$loglik
      }

      if (q > 0 && table[p + 1, q + 1] + eps_tol < table[p + 1, q]) {
        is_consistent = FALSE
        break
      } else if (p > 0 && table[p + 1, q + 1] + eps_tol < table[p, q + 1]) {
        is_consistent = FALSE
        break
      }

    }

    if (!is_consistent) break
  }

  is_consistent
}

