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
#' @param max_iters Maximum number of random restarts for methods "CSS-ML" and
#'    "ML". If set to 1, the results of this algorithm is the same as
#'    [stats::arima()] if argument \code{diffuseControl} is also set as TRUE.
#'    \code{max_iters} is often not reached because the condition
#'    \code{max_repeats} is typically achieved first.
#' @param max_repeats Integer. If the last \code{max_repeats} random starts did
#'    not result in improved likelihoods, then stop the search. Each result of
#'    the optim function is only considered to improve the likelihood if it does
#'    so by more than \code{eps_tol}.
#' @param eps_tol Tolerance for accepting a new solution to be better than a
#'    previous solution. The default corresponds to a one ten-thousandth
#'    unit increase in log-likelihood.
#'
#' @return A matrix containing the model AIC values.
#' @export
#' @examples
#' aicTable(presidents, 3, 2)
aicTable <- function(data, P, Q, D = 0, max_repeats = 10, max_iters = 100, eps_tol = 1e-4){

  if (!is.numeric(P) | !is.numeric(Q) | !is.numeric(D)) {
    stop("'P', 'Q' and 'D' must be numeric.")
  }

  P <- as.integer(P)
  Q <- as.integer(Q)
  D <- as.integer(D)

  table <- matrix(NA, (P + 1), (Q + 1))
  for(p in 0:P) {
    for(q in 0:Q) {
      table[p + 1, q + 1] <- arima(data, order = c(p, D, q), max_iters = max_iters)$aic
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
#' @param data a time series object, or a dataset that can be used as input into
#'    the [arima] function.
#' @param P a positive integer value representing the maximum number of AR
#'    coefficients that should be included in the table.
#' @param Q a positive integer value representing the maximum number of MA
#'    coefficients that should be included in the table.
#' @param D a positive integer value representing the degree of differencing
#' @param max_iters Maximum number of random restarts for methods "CSS-ML" and
#'    "ML". If set to 1, the results of this algorithm is the same as
#'    [stats::arima()] if argument \code{diffuseControl} is also set as TRUE.
#'    \code{max_iters} is often not reached because the condition
#'    \code{max_repeats} is typically achieved first.
#' @param max_repeats Integer. If the last \code{max_repeats} random starts did
#'    not result in improved likelihoods, then stop the search. Each result of
#'    the optim function is only considered to improve the likelihood if it does
#'    so by more than \code{eps_tol}.
#' @param eps_tol Tolerance for accepting a new solution to be better than a
#'    previous solution. The default corresponds to a one ten-thousandth
#'    unit increase in log-likelihood.
#' @param method string that must be "arima2" or "stats", indicating which
#'    package should be used to fit the ARIMA model.
#'
#' @return Boolean. True if the table is consistent in the sense that larger
#'    models have likelihood that is greater than or equal to all smaller
#'    models.
#' @noRd
#'
#' @examples .checkTable(presidents, 3, 2)
.checkTable <- function(data, P, Q, D = 0, max_repeats = 10, max_iters = 100, eps_tol = 1e-4,
                        method = 'arima2') {

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
        table[p + 1, q + 1] <- arima(data, order = c(p, D, q), max_iters = max_iters, max_repeats = max_repeats, eps_tol = eps_tol)$loglik
      } else {
        table[p + 1, q + 1] <- stats::arima(data, order = c(p, D, q))$loglik
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

