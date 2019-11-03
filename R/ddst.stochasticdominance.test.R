#' Data Driven Two-Sample Test Against Stochastic Dominance
#'
#' Performs data driven smooth non-parametric two-sample test against one-sided alternatives
#' (stochastic dominance).
#' Suppose that we have random samples from two distributions F and G.
#' The null hypothesis is that F(x) < G(x) for some x while the alternative is that at
#' F(x) >= G(x) for all x with strict inequality for at least one x.
#' Detailed description of the test statistic is provided in Ledwina and Wylupek (2012).
#'
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param t a positive number, penalty for model selection rule, see package description
#' @param d an integer, number of coordinates that measure potential deviation from null hypothesis
#' @param alpha a positive number, significance level for test statistic. Note that tuning parameter t depends on alpha. See package description
#' @param ... further arguments
#'
#' @references Two-sample test against one-sided alternatives. Ledwina and Wylupek (2012). \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9469.2011.00787.x}
#' @export
#' @examples
#' # H0 is true
#' x <- runif(80)
#' y <- runif(80)
#' t <- ddst.againststochdom.test(x, y)
#' t
#' plot(t)
#'
#' # known fixed alternative
#' x <- runif(80)
#' y <- rbeta(80,4,2)
#' t <- ddst.againststochdom.test(x, y)
#' t
#' plot(t)
#' @keywords htest
`ddst.againststochdom.test` <-
  function(x,
           y,
           t = 2.2,
           d = 4,
           alpha = 0.05,
           ...) {
    coord = ddst.stochasticdominance.Nk(x, y, t = t, k.N = d, alpha = alpha)    # coord square times n

    l = coord$T
    attr(l, "names") = "T"
    vt = coord$V.T
    attr(vt, "names") = "VT"
    result = list(statistic = vt,
                  parameter = l,
                  coordinates = coord$L,
                  method = "Data Driven Test Against Stochastic Dominance")
    result$data.name = paste(paste(as.character(substitute(x)), collapse = ""),
                             ", alpha: ",
                             alpha,
                             ", t: ",
                             t,
                             ", Dmax: ",
                             d,
                             sep = "")
    class(result) = c("htest", "ddst.test", "ddst.stochasticdominance.test")

    result
  }
