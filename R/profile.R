#' Profile for `Arima2` object
#'
#' This function performs profile log-likelihood of an `Arima2` function.
#'
#' The parameter `which` specifies parameter in the following vector will be
#' profiled over:
#' \deqn{\phi_1, \ldots, \phi_p, \theta_1, \ldots, \theta_q, \Phi_1, \ldots, \Phi_P, \Theta_1, \ldots, \Theta_Q, \mu}
#' where \eqn{p, q} are non-negative integers representing the number of AR and
#' MA coefficients of `fitted`, respectively, and \eqn{\phi_i} are the AR
#' coefficients, \eqn{\theta_i} are the MA coefficients, \eqn{\Phi_i} are the
#' seasonal AR coefficients, \eqn{\Theta_i} are the seasonal MA coefficients and
#' \eqn{\mu} is the model intercept.
#'
#' @param fitted An `Arima2` object that has been fit to data.
#' @param d Integer number of differences. Should match the differences used to
#'    obtain the `fitted` object.
#' @param npts Integer number of points to evaluate the profile.
#' @param lower Numeric lower bound for the profile search.
#' @param upper Numeric upper bound for the profile search.
#' @param which Integer indicating which parameter to perform the profile over.
#'    See Details section for more information.
#' @param ... additional arguments needed for the profile function
#'
#'
#' @importFrom methods hasArg
#' @return data.frame object containing the results of the profile likelihood.
#' @export
profile.Arima2 <- function(fitted, d = 0, npts = 100L, lower = -1, upper = 1, which = 1L, ...) {

  if (!is.numeric(which)) stop("argument `which` must be numeric.")
  if (!is.numeric(lower) | !is.numeric(upper)) stop('arguments `lower` and `upper` must both be numeric.')
  if (lower >= upper) stop("argument `lower` must be less than `upper`.")
  if (!is.numeric(npts)) stop("argument `npts` must be numeric.")

  coefs <- fitted[['coef']]

  has_seasonal <- FALSE
  if (any(grepl("sar\\d+|sma\\d+", names(coefs)))) {
    has_seasonal <- TRUE
    if (!hasArg(period)) {
      period <- 1
      warning("Fitted model has seasonal component but no `period` was given. Setting `period` to 1.")
    }

    if (!hasArg(seas_d)) {
      seas_d <- 0
      warning("Fitted model has seasonal component but no `seas_d` was given. Setting `seas_d` to 0.")
    }
  }

  which <- as.integer(which)
  if (which > length(coefs)) stop("argument `which` must be less than the number of parameters in `fitted`.")

  npts <- as.integer(npts)

  pll_df <- matrix(nrow = npts, ncol = length(coefs) + 1L)
  colnames(pll_df) <- c(names(coefs), "loglik")
  pll_df[, which] <- seq(from = lower, to = upper, length.out = npts)

  fitted_order <- c(
    sum(grepl("^ar\\d+$", names(coefs))),
    d,
    sum(grepl("^ma\\d+$", names(coefs)))
  )

  if (has_seasonal) {
    fitted_seasonal <- list(
      order = c(
        sum(grepl("^sar\\d+$", names(coefs))),
        seas_d,
        sum(grepl("^sma\\d+$", names(coefs)))
      ),
      period = period
    )
  }

  for (i in 1:npts) {
    tmp_fixed <- pll_df[i, -(length(coefs) + 1)]

    if (has_seasonal) {
      tmp_out <- tryCatch(arima(
        fitted$x, order = fitted_order, seasonal = fitted_seasonal,
        fixed = tmp_fixed, transform.pars = FALSE, max_iters = 1),
        error = function(e) list(coef = tmp_fixed, loglik = NA))
    } else {
      tmp_out <- tryCatch(arima(
        fitted$x, order = fitted_order,
        max_iters = 1, fixed = tmp_fixed, transform.pars = FALSE
      ),
      error = function(e) list(coef = tmp_fixed, loglik = NA))
    }

    pll_df[i, ] <- c(tmp_out$coef, tmp_out$loglik)
  }

  as.data.frame(pll_df)
}
