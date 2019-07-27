#' Data Driven k-Sample Test
#'
#' Performs data driven smooth test for the classical two- and k-sample problems as described in details in Wylupek (2010).
#'
#' @param x a list with vectors or vectors (in thos case other vectors will be in ...)
#' @param d_N an integer, number of coordinates, see package description
#' @param c a positive number, penalty for model selection rule, see package description
#' @param B an integer specifying the number of replicates used in p-value computation
#' @param compute.p a logical value indicating whether to compute a p-value
#' @param ... other vectors
#'
#' @export
#'
#' @examples
#' # H0 is true
#' x = runif(80)
#' y = runif(80)
#' z = runif(80)
#' ddst.ksample.test(x, y, z, compute.p=TRUE)
#' ddst.ksample.test(list(x, y, z), compute.p=TRUE)
#'
#' # known fixed alternative
#' x = runif(80)
#' y = rbeta(80,4,2)
#' w = rnorm(30)
#' z = rexp(10, 1)
#' ddst.ksample.test(x, y, w, z, compute.p=TRUE)
#' ddst.ksample.test(list(x, y, w, z), compute.p=TRUE)
#' @keywords htest
`ddst.ksample.test` <-
  function(x,
           ...,
           d_N = 12, c = 2.3,
           B = 1000,
           compute.p = FALSE) {
    if (is.list(x)) {
      # x is list with coordinates
      x.vector <- unlist(x)
      n <- sapply(x, length)
    } else {
      # x is a vector, other vectors are in ...
      x.vector <- unlist(c(x, list(...)))
      n <- c(length(x),
             sapply(list(...), length))
    }
    coord = ddst.ksample.Nk(x.vector, n, d_N = d_N, c = c)

    l = coord[1]
    attr(l, "names") = "W.T"
    t = NA
    attr(t, "names") = "l"
    result = list(statistic = l,
                  parameter = t,
                  method = "Data Driven k-Sample Test")
    class(result) = "htest"

    if (compute.p) {
      tmp = numeric(B)
      for (i in 1:B) {
        xy <- sample(x.vector)

        tmp[i] = ddst.ksample.Nk(xy, n, d_N = d_N, c = c)
      }
      result$p.value = mean(tmp > l)
    }
    result
  }
