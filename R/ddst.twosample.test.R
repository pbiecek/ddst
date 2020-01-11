#' Data Driven Smooth Test for Two-Sample Problem
#'
#' Performs data driven smooth test for the classical two-sample problem.
#' It is  aspecial case of data driven test for k-samples.
#' Detailed description of the test statistic is provided in Wylupek (2010).
#'
#' @param x a (non-empty) numeric vector of data
#' @param y a (non-empty) numeric vector of data
#' @param d.n an integer specifying the maximum dimension considered, only for advanced users
#' @param c a calibrating parameter in the penalty in the model selection rule
#' @param B an integer specifying the number of runs for a p-value and a critical value computation if any
#' @param compute.p a logical value indicating whether to compute a p-value or not
#' @param alpha a significance level
#' @param compute.cv a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not
#'
#' @references Data-driven k-sample tests. Wylupek (2010) \url{https://www.jstor.org/stable/40586684?seq=1}
#' @export
#' @examples
#' set.seed(7)
#' # H0 is true
#' x <- runif(80)
#' y <- runif(80)
#' t <- ddst.twosample.test(x, y)
#' t
#' plot(t)
#'
#' # H0 is false
#' x <- runif(80)
#' y <- rexp(80, 1)
#' t <- ddst.twosample.test(x, y)
#' t
#' plot(t)
#'
#' @keywords htest
`ddst.twosample.test` <-
  function(x,
           y,
           d.n = 12,
           c = 2,
           B = 10000,
           compute.p = FALSE,
           alpha = 0.05,
           compute.cv = FALSE) {
    res <- ddst.ksample.test(list(x, y), d.n = d.n, c = c, B = B, compute.p = compute.p,
                             alpha = alpha, compute.cv = compute.cv)
    res$coordinates <- res$coordinates[1,]
    res$method <- "Data Driven Smooth Test for Two-Sample Problem"
    class(res) = c("htest", "ddst.test", "ddst.twosample.test")
    res
  }
