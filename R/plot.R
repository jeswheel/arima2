#' Plot `Arima2` object
#'
#' @param x An `Arima2` object.
#' @param title Title of plot
#' @param tick.lab Time vector
#' @param ... Other parameters
#'
#' @return plot
#' @export plot.arima2
#'
#' @examples
#' mod <- plot.arima2(arima2(lh, order = c(1,0,1)))
plot.arima2 <- function(x, title = NULL, tick.lab = NULL, ...) {

  ### check if tick.lab is left as default
  if(is.null(tick.lab)){
    # x is an ARIMA2 object, x$x contains the data.
    xaxis <- 1:length(x$x)
    yaxis <- x$x
    df <- data.frame(xaxis = xaxis, yaxis = yaxis) ## ggplot2 needs a data frame of the data
    plot <- ggplot2::ggplot(df, ggplot2::aes(x = xaxis, y = yaxis)) +
      ggplot2::geom_line() +
      ggplot2::scale_x_continuous(breaks = seq(min(df$x), max(df$x), by = 1))

    if(!is.null(title)){
      plot + ggplot2::labs(title = title)
      }
  }


  ### check if tick.lab contains numeric values (could be years)
  else if(is.numeric(tick.lab)){
    ## important to add because current messages cause debugging difficulty
    if(length(x$x) != length(tick.lab))stop("The length of data points does not equal the length of time points. Check \"tick.lab\" parameter.")
    # x is an ARIMA2 object, x$x contains the data.
    xaxis <- tick.lab
    yaxis <- x$x
    df <- data.frame(xaxis = xaxis, yaxis = yaxis) ## ggplot2 needs a data frame of the data
    plot <- ggplot2::ggplot(df, ggplot2::aes(x = xaxis, y = yaxis)) +
      ggplot2::geom_line() +
      ggplot2::scale_x_continuous(breaks = seq(min(df$x), max(df$x), by = 1))

    if(!is.null(title)){
      plot + ggplot2::labs(title = title)
    }
  }

  ### check if tick.lab contains characters
  else if(!is.numeric(tick.lab)){
    ## important to add because current messages cause debugging difficulty
    if(length(x$x) != length(tick.lab))stop("The length of data points does not equal the length of time points. Check \"tick.lab\" parameter.")
    # x is an ARIMA2 object, x$x contains the data.
    xaxis <- tick.lab
    yaxis <- x$x
    df <- data.frame(xaxis = xaxis, yaxis = yaxis) ## ggplot2 needs a data frame of the data

    plot <- ggplot2::ggplot(df, ggplot2::aes(x = xaxis, y = yaxis)) +
      ggplot2::geom_line(ggplot2::aes(group = 1))

    if(!is.null(title)){
      plot + ggplot2::labs(title = title)
    }
  }


  return(plot)
}

