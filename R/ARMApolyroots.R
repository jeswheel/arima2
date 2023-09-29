#' ARMA polyroots
#'
#' This function calculates the roots of the AR or MA polynomials that
#' correspond to an ARMA model.
#'
#' @param model Either of fitted object of class `Arima` (i.e., the output of
#'     either [stats::arima()] or [arima]), a list with named elements at
#'     least one of the named elements `ar` or `ma`, or a vector with named
#'     elements, such as `c("ar1" = 0.3, "ar2" = -0.2, "ma1" = 0.14)` Seasonal
#'     coefficients are ignored.
#' @param type character of value "AR" or "MA", indicating whether or not the
#'     AR or MA polynomial roots are desired.
#'
#' @returns A numeric vector containing the roots of the MA or AR polynomials
#' @export
#'
#' @examples
#' set.seed(123456)
#' ARMApolyroots(sample_ARMA_coef((order = c(2, 1))), type = "AR")
#'
#' mod <- arima(lh, order = c(3,0,0))
#' ARMApolyroots(mod, type = "AR")
ARMApolyroots <- function(model, type = c("AR", "MA")) {
  type <- match.arg(type)

  if (methods::is(model, "numeric")) {
    model <- list(
      ar = model[grepl("^ar[[:digit:]]+$", names(model))],
      ma = model[grepl("^ma[[:digit:]]+$", names(model))]
    )
  } else if (inherits(model, "Arima")) {
    model <- list(
      ar = model$coef[grepl("^ar[[:digit:]]+$", names(model$coef))],
      ma = model$coef[grepl("^ma[[:digit:]]+$", names(model$coef))]
    )
  }

  if (type == "AR") {
    if (length(model$ar) == 0L) {
      stop("Roots of AR polynomial requested, but no AR coefficients provided")
    }

    p <- max(which(c(1, -model$ar) != 0)) - 1
    return(polyroot(c(1, -model$ar[1L:p])))

  } else if (type == "MA") {
    if (length(model$ma) == 0L) {
      stop("Roots of MA polynomial requested, but no MA coefficients provided")
    }

    p <- max(which(c(1, model$ma) != 0)) - 1
    return(polyroot(c(1, model$ma[1L:p])))
  }
}

