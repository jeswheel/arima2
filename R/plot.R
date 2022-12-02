#' Plot `Arima2` object
#'
#' This function plots an `Arima2` object.
#'
#' TODO: This could plot data, data + sims, or roots, roots + data, roots + data + sims.
#'
#' @param object An `Arima2` object.
#'
#' @return NULL
#' @export
#'
#' @examples
#' mod <- plot(arima2(lh, order = c(1,0,1)))
plot.Arima2 <- function(object) {
  plot(object$x, type = 'l')
}
