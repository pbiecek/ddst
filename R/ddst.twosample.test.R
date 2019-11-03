#' Data Driven Two-Sample Test
#'
#' Performs data driven smooth test for the classical two-sample problem.
#' It is  aspecial case of data driven test for k-samples.
#' Detailed description of the test statistic is provided in Wylupek (2010).
#'
#' @param x a list with two vectors or a single vector.
#' @param ... if x is a single vector, then the second vector is provided in the ...
#' @param d_N an integer, number of coordinates that measure potential deviation from null hypothesis
#' @param c a positive number, penalty for model selection rule. Section 5.1 in Wylupek (2010), suggest that a good choice is c = 2, when k = 2, and  c = 2.3, when k >= 3.
#'
#' @references Data-driven k-sample tests. Wylupek (2010) \url{https://www.jstor.org/stable/40586684?seq=1}
#' @export
#' @examples
#' set.seed(7)
#' # H0 is true
#' x <- runif(80)
#' y <- runif(80)
#' t <- ddst.twosample.test(x, y)
#' t <- ddst.twosample.test(list(x, y))
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
           ...,
           d_N = 12, c = 2) {
    res <- ddst.ksample.test(x, ..., d_N = d_N, c = c)
    res$coordinates <- res$coordinates[1,]
    res$method <- "Data Driven Two-Sample Test"
    class(res) = c("htest", "ddst.test", "ddst.twosample.test")
    res
  }
