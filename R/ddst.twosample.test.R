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
#' z = runif(80)
#' ddst.uniform.test(z, compute.p=TRUE)
#'
#' # known fixed alternative
#' z = rnorm(80,10,16)
#' ddst.uniform.test(pnorm(z, 10, 16), compute.p=TRUE)
#'
#' # H0 is false
#' z = rbeta(80,4,2)
#' (t = ddst.uniform.test(z, compute.p=TRUE))
#' t$p.value
#' @keywords htest
`ddst.twosample.test` <-
  function(x,
           y,
           t = 2.2,
           k.N = 4,
           B = 1000,
           compute.p = FALSE,
           ...) {
#    n = length(x)
#    if (n < 5)
#      stop("length(x) should be at least 5")


    coord = ddst.twosample.Nk(x, y, t = t, k.N = k.N)    # coord square times n

    l = ddst.IIC(coord, n, c)
    attr(l, "names") = "n. coord"
    t = coord[l]
    attr(t, "names") = "WT"
    result = list(statistic = t,
                  parameter = l,
                  method = "Data Driven Smooth Test for Uniformity")
    result$data.name = paste(paste(as.character(substitute(x)), collapse =
                                     ""),
                             ",   base: ",
                             method.name,
                             "   c: ",
                             c,
                             sep = "")
    class(result) = "htest"
    if (compute.p) {
      tmp = numeric(B)
      for (i in 1:B) {
        y = runif(n)
        tmpC = ddst.uniform.Nk(y, base, Dmax = Dmax)
        l = ddst.IIC(tmpC, n, c)
        tmp[i] = tmpC[l]
      }
      result$p.value = mean(tmp > t)
    }
    result
  }
