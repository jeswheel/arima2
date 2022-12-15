#' Plot `Arima2` object
#'
#' This function plots time series data loaded from an `Arima2` object or plots inverse roots of the AR or MA polynomials in a fitted ARIMA model on the complex unit circle. Note: Layers may be added to `Arima2` plots using `ggplot2`'s plus operator.
#'
#' @param x An `Arima2` object. This parameter is an object created using the function `arima2()`.
#' @param title Title of plot
#' @param tick.lab Time vector of numeric or character/string type.
#' @param roots Would you instead prefer to plot the roots on a unit circle? Insert logical type here.
#' @param ... Other parameters
#'
#' @return `Arima 2` plot. Type of plot is indicated through `roots` parameter.
#' @export plot.arima2
#'
#' @examples
#' mod <- plot.arima2(arima2(lh, order = c(1,0,1)))
#' my.plot <- plot.arima2(x = arima2(lh, order = c(3,0,1)), roots = TRUE)
plot.arima2 <- function(x, roots = FALSE, title = NULL, tick.lab = NULL, ...) {

  if(!is.null(title) & !is.character(title))stop("'title' should be character type.")
  if(!is.logical(roots))stop("'roots' should be logical type.")
  if(roots == FALSE){
    ### check if tick.lab is left as default
    if(is.null(tick.lab)){
      # x is an ARIMA2 object, x$x contains the data.
      xaxis <- 1:length(x$x)
      yaxis <- x$x
      df <- data.frame(xaxis = xaxis, yaxis = yaxis) ## ggplot2 needs a data frame of the data
      plot <- ggplot2::ggplot(df, ggplot2::aes(x = xaxis, y = yaxis)) +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = seq(min(df$x), max(df$x), by = 1)) +
        ggplot2::labs(x = "Time", y = "Value")

      if(!is.null(title)){
        plot$labels$title = title
        }
    }


    ### check if tick.lab contains numeric values (could be years)
    else if(is.numeric(tick.lab)){
      ## important to add because current messages cause debugging difficulty
      if(length(x$x) != length(tick.lab))stop("The length of data points does not equal the length of time points. Check 'tick.lab' parameter.")
      # x is an ARIMA2 object, x$x contains the data.
      xaxis <- tick.lab
      yaxis <- x$x
      df <- data.frame(xaxis = xaxis, yaxis = yaxis) ## ggplot2 needs a data frame of the data
      plot <- ggplot2::ggplot(df, ggplot2::aes(x = xaxis, y = yaxis)) +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(breaks = seq(min(df$x), max(df$x), by = 1)) +
        ggplot2::labs(x = "Time", y = "Value")

      if(!is.null(title)){
        plot$labels$title = title
      }
    }

    ### check if tick.lab contains characters
    else if(!is.numeric(tick.lab)){
      ## important to add because current messages cause debugging difficulty
      if(length(x$x) != length(tick.lab))stop("The length of data points does not equal the length of time points. Check 'tick.lab' parameter.")
      if(!is.numeric(tick.lab) & !is.character(tick.lab))stop("'tick.lab' must be a numeric or character type.")
      # x is an ARIMA2 object, x$x contains the data.
      xaxis <- as.factor(tick.lab)
      levels(xaxis) <- tick.lab ## prevent reordering xaxis in ggplot
      yaxis <- x$x
      df <- data.frame(xaxis = xaxis, yaxis = yaxis) ## ggplot2 needs a data frame of the data

      plot <- ggplot2::ggplot(df, ggplot2::aes(x = as.factor(xaxis), y = yaxis)) +
        ggplot2::geom_line(ggplot2::aes(group = 1)) +
        ggplot2::labs(x = "Time", y = "Value")

      if(!is.null(title)){
        plot$labels$title = title
      }
    }
  }

  ## displays all inverse roots within unit circle
  else{
    if(!is.null(title) | !is.null(tick.lab))warning("Parameters 'title' and 'tick.lab' will have no effect. See 'ggplot2' package documentation for more information.")
    plot <- forecast::autoplot(x) +
      ggplot2::labs(color = "Unit Circle") +
      ggplot2::theme(axis.title = ggplot2::element_text(face="bold"))
    plot$layers[[4]]$aes_params$size = 1.5 ## making points smaller
    }


  return(plot)
}

