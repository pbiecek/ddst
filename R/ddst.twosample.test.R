#' Data Driven Two-Sample Test Against One-Sided Alternatives
#'
#' Performs data driven smooth non-parametric two-sample test against one-sided alternatives (stochastic dominance)
#'
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values
#' @param t a positive number, penalty for model selection rule, see package description
#' @param k.N (TODO: is it D ?)
#' @param B an integer specifying the number of replicates used in p-value computation
#' @param compute.p a logical value indicating whether to compute a p-value
#' @param ... further arguments
#'
#' @export
#'
#' @examples
#' # H0 is true
#' x = runif(80)
#' y = runif(80)
#' ddst.twosample.test(x, y, compute.p=TRUE)
#'
#' # known fixed alternative
#' x = runif(80)
#' y = rbeta(80,4,2)
#' ddst.twosample.test(x, y, compute.p=TRUE)
#' @keywords htest
`ddst.twosample.test` <-
  function(x,
           y,
           t = 2.2,
           k.N = 4,
           B = 1000,
           compute.p = FALSE,
           alpha = 0.05,
           ...) {
    coord = ddst.twosample.Nk(x, y, t = t, k.N = k.N, alpha = alpha)    # coord square times n

    l = coord[3]
    attr(l, "names") = "T"
    t = coord[1]
    attr(t, "names") = "V.T"
    result = list(statistic = t,
                  parameter = l,
                  method = "Data Driven Two-Sample Test")
    class(result) = "htest"

    if (compute.p) {
      tmp = numeric(B)
      for (i in 1:B) {
        xy <- sample(c(x,y))
        new.x <- xy[1:length(x)]
        new.y <- xy[(length(x)+1):length(xy)]

        tmp[i] = ddst.twosample.Nk(x, y, t = t, k.N = k.N, alpha = alpha)[1]
      }
      result$p.value = mean(tmp > t)
    }
    result
  }
