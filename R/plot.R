#' Plot `Arima2` object
#'
#' This function plots time series data loaded from an `Arima2` object or plots
#' inverse roots of the AR or MA polynomials in a fitted ARIMA model on the
#' complex unit circle.
#'
#' The output of this function is a `ggplot` object, so layers may be added to
#' the output of this function using `ggplot2`'s plus operator.
#'
#' @param x An `Arima2` object. This parameter is an object created using the function `arima2()`.
#' @param title Title of plot
#' @param tick.lab Time vector of numeric or character/string type.
#' @param roots Would you instead prefer to plot the roots on a unit circle? Insert logical type here.
#' @param ... Other parameters
#'
#' @import ggplot2
#' @return `Arima 2` plot, which is a `ggplot2` object. Type of plot is
#'    indicated through `roots` parameter.
#' @export
#'
#' @examples
#' plot(arima(lh, order = c(1,0,1)))
#' plot(x = arima(lh, order = c(3,0,1)), roots = FALSE)
plot.Arima2 <- function(x, roots = TRUE, title = NULL, tick.lab = NULL, ...) {

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

    ar_inv_roots <- 1 / ARMApolyroots(x)
    ma_inv_roots <- 1 / ARMApolyroots(x, type = "MA")

    ar_df <- data.frame(
      root_type = "AR",
      value = ar_inv_roots
    )

    ma_df <- data.frame(
      root_type = "MA",
      value = ma_inv_roots
    )

    all_roots <- rbind(ar_df, ma_df)
    all_roots$Real <- Re(all_roots$value)
    all_roots$Imaginary <- Im(all_roots$value)
    all_roots$inside <- ifelse(Mod(all_roots$value) < 1, "inside", "outside")

    if(!is.null(title) | !is.null(tick.lab))warning("Parameters 'title' and 'tick.lab' will have no effect. See 'ggplot2' package documentation for more information.")

    plot <- ggplot(all_roots, aes(x = .data[["Real"]], y = .data[["Imaginary"]], col = .data[["inside"]])) +
      annotate("path", x = cos(seq(0, 2*pi, length.out = 100)),
               y = sin(seq(0, 2*pi, length.out = 100))) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point() +
      facet_wrap(
        ~root_type,
        nrow = 1,
        labeller = as_labeller(
          c("AR" = "Inverse AR roots", "MA" = "Inverse MA roots"))
      ) +
      coord_fixed() +
      labs(color = "Unit Circle") +
      theme_bw() +
      theme(axis.title = ggplot2::element_text(face="bold"),
            legend.position = 'bottom')

    plot$layers[[4]]$aes_params$size = 1.5 ## making points smaller
  }


  return(plot)
}

