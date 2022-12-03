#' Plot `Arima2` object
#'
#' This function plots an `Arima2` object.
#'
#' TODO: This could plot data, data + sims, or roots, roots + data, roots + data + sims.
#'
#' @param x An `Arima2` object.
#' @param ... All inputs to the plot function
#' @return NULL
#' @export
#'
#' @examples
#' mod <- plot(arima2(lh, order = c(1,0,1)))
plot.Arima2 <- function(x, ...) {
  # x is an ARIMA2 object, x$x contains the data.
  xaxis <- 1:length(x$x)
  yaxis <- x$x
  plot(xaxis, yaxis, type = 'l', ...)
}
