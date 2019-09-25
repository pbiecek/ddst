#' Plot Function fo Data Driven Tests
#'
#' Plots coordinates for selected test statistic
#'
#' @param x result from ddst test function
#'
#' @export
#' @examples
#' # H0 is true
#' x <- runif(80)
#' y <- runif(80)
#' z <- runif(80)
#' t <- ddst.ksample.test(x, y, z)
#' plot(t)
`plot.ddst.test` <-
  function(x) {
    coordinates <- x$coordinates
    if (is.null(dim(coordinates))) {
      # coordinates is a vector
      # coordinates is a matrix
      val <- c(coordinates)
      pos <- 1:length(coordinates)
      df <- data.frame(val, pos)
      pl <- ggplot(df, aes(factor(pos), val)) +
        geom_col() +
        xlab("coordinate") + ylab("")
    } else {
      # coordinates is a matrix
      val <- c(coordinates)
      vec <- rep(1:nrow(coordinates), each = ncol(coordinates))
      pos <- rep(1:ncol(coordinates), times = nrow(coordinates))
      df <- data.frame(val, vec, pos)
      pl <- ggplot(df, aes(factor(pos), val)) +
        geom_col() + facet_grid(vec~.) +
        xlab("coordinate") + ylab("")
    }
    pl +
      theme_bw()
  }
