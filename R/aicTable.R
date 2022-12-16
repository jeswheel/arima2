#' ARIMA AIC table
#'
#'Construct table of AIC for all combinations 0<=p<=P and 0<=q<=Q
#'
#' This function creates an AIC table for ARMA models of varying sizes.
#' Each row for the table corresponds to a different AR value, and each column
#' of the table corresponds to a different MA value.
#'
#' @param data a time series object, or a dataset that can be used as input into
#'    the [arima2] function.
#' @param P a positive integer value representing the maximum number of AR
#'    coefficients that should be included in the table.
#' @param Q a positive integer value representing the maximum number of MA
#'    coefficients that should be included in the table.
#' @param D a positive integer value representing the degree of differencing
#'    to be used for each \eqn{ARMA(p, q)} model.
#' @param nrestart number of restarts for the [arima2] function.
#'
#' @export
#' @examples
#' aicTable(presidents, 3, 2)
aicTable <- function(data, P, Q, D = 0, nrestart = 10){

  if (!is.numeric(P) | !is.numeric(Q) | !is.numeric(D)) {
    stop("'P', 'Q' and 'D' must be numeric.")
  }

  P <- as.integer(P)
  Q <- as.integer(Q)
  D <- as.integer(D)

  table <- matrix(NA, (P + 1), (Q + 1))
  for(p in 0:P) {
    for(q in 0:Q) {
      table[p + 1, q + 1] <- arima2(data, order = c(p, D, q), nrestart = nrestart)$aic
    }
  }
  dimnames(table) <- list(paste("AR", 0:P, sep = ""), paste("MA", 0:Q, sep = ""))
  table
}
