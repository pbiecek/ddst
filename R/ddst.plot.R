#' Plot Function fo Data Driven Tests
#'
#' Plots coordinates for selected test statistics
#'
#' @param x result from ddst test function
#' @param ... currently not used
#'
#' @importFrom stats pexp pnorm qnorm rexp rnorm runif
#' @importFrom ggplot2 aes facet_grid geom_col ggplot ggtitle theme_bw xlab ylab
#' @export
#' @examples
#' # H0 is true
#' x <- runif(80)
#' y <- runif(80)
#' z <- runif(80)
#' t <- ddst.ksample.test(x, y, z)
#' plot(t)
`plot.ddst.test` <-
  function(x, ...) {
    coordinates <- x$coordinates
    pos <- val <- vec <- NULL # keep ggplot2 clean

    if (is.null(dim(coordinates))) {
      # coordinates is a vector
      # coordinates is a matrix
      val <- c(coordinates)
      if (is.null(names(coordinates))) {
        pos <- 1:length(coordinates)
      } else {
        pos <- names(coordinates)
      }
      df <- data.frame(val, pos)
      pl <- ggplot(df, aes(factor(pos), val)) +
        geom_col() +
        xlab("j") + ylab("")
    } else {
      # coordinates is a matrix
      val <- c(coordinates)
      vec <- rep(1:nrow(coordinates), each = ncol(coordinates))
      if (is.null(colnames(coordinates))) {
        pos <- rep(1:ncol(coordinates), times = nrow(coordinates))
      } else {
        pos <- rep(colnames(coordinates), times = nrow(coordinates))
      }
      df <- data.frame(val, vec, pos)
      pl <- ggplot(df, aes(factor(pos), val)) +
        geom_col() + facet_grid(vec~.) +
        xlab("j") + ylab("")
    }
    pl +
      theme_bw() +
      ggtitle(paste0("Method: ", x$method),
              paste0( "Test statistic ", attr(x$statistic, "names"), ": ",
                      signif(x$statistic, 3)))
  }
