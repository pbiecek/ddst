#' Data Driven Test for k-Sample Umbrella Alternatives
#'
#' Performs data driven smooth test for so-called umbrella alternatives in k-sample problem;
#' see Wylupek (2016).
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
#' ddst.umbrella.test(x, y, z, compute.p=TRUE)
#' ddst.umbrella.test(x, y, z, compute.p=TRUE, type = "S")
#' ddst.umbrella.test(x, y, z, compute.p=TRUE, type = "M")
#' ddst.umbrella.test(list(x, y, z), compute.p=TRUE)
#'
#' # known fixed alternative
#' x = runif(80)
#' y = rbeta(80,4,2)
#' w = rnorm(30)
#' z = rexp(10, 1)
#' ddst.umbrella.test(x, y, w, z)
#' ddst.umbrella.test(x, y, w, z, type = "S")
#' ddst.umbrella.test(x, y, w, z, type = "M")
#' ddst.umbrella.test(list(x, y, w, z))
#' @keywords htest
`ddst.umbrella.test` <-
  function(x,
           ...,
           tlh.p = 2.2, tl.n = 2.2, p = 3,
           type = "T",
           B = 100,
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

    # type must be S, T or M
    index <- c(T = 2, S = 1, M = 3)[type]
    coord = ddst.umbrella.Nk(x.vector, n, tlh.p = tlh.p, tl.n = tl.n, p = p)
    # TODO: S, T, M - which test to use?
    # what about default for tlh.p
    # now all three are being calculated

    l = coord[index]
    attr(l, "names") = type
    t = NA
    attr(t, "names") = "l"
    result = list(statistic = l,
                  parameter = t,
                  method = "Data Driven k-Sample Umbrella Test")
    class(result) = "htest"

    if (compute.p) {
      tmp = numeric(B)
      for (i in 1:B) {
        xy <- sample(x.vector)

        tmp[i] = ddst.umbrella.Nk(xy, n, tlh.p = tlh.p, tl.n = tl.n, p = p)[index]
      }
      result$p.value = mean(tmp > l)
    }
    result
  }
