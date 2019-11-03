#' Data Driven Test for k-Sample Umbrella Alternatives
#'
#' Performs data driven smooth test for so-called umbrella alternatives in k-sample problem.
#' Suppose that we have random samples from k distributions F_i where i = 1, ..., k.
#' The null hypothesis is that there is no umbrella pattern,
#' i.e. F_1 >= ... >= F_p <= ... <= F_k and F_i != F_j for some i and j.
#' The alternative is that there is an umbrella pattern
#' i.e. F_1 >= ... >= F_p <= ... <= F_k and F_i != F_j for some i and j.
#' Detailed description of the test statistic is provided in Wylupek (2016).
#'
#' @param x a list with k-vectors or a single vector.
#' @param ... if x is a single vector, then remaing k-1 vectors are provided in the ... argument
#' @param tlh.p a positive number, penalty for model selection rule
#' @param tl.n a positive number, penalty for model selection rule
#' @param p peak of the umbrella
#'
#' @references An automatic test for the umbrella alternatives. Wylupek (2016) \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/sjos.12231}
#' @export
#' @examples
#' set.seed(7)
#' # H0 is true
#' x = runif(80)
#' y = runif(80) + 0.2
#' z = runif(80)
#' ddst.umbrella.test(x, y, z, p = 2)
#' t <- ddst.umbrella.test(list(x, y, z), p = 2)
#' t
#' plot(t)
#'
#' # known fixed alternative
#' x1 = rnorm(80)
#' x2 = rnorm(80) + 2
#' x3 = rnorm(80) + 4
#' x4 = rnorm(80) + 3
#' x5 = rnorm(80) + 2
#' x6 = rnorm(80) + 1
#' x7 = rnorm(80)
#' t <- ddst.umbrella.test(list(x1, x2, x3, x4, x5, x6, x7), p = 3)
#' t
#' plot(t)
#'
#' t <- ddst.umbrella.test(list(x1, x2, x3, x4, x5, x6, x7), p = 5)
#' t
#' plot(t)
#'
#' @keywords htest
`ddst.umbrella.test` <-
  function(x,
           ...,
           tlh.p = 2.2, tl.n = 2.2,
           p = 3) {
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

    coord = ddst.umbrella.Nk(x.vector, n,
                             tlh.p = tlh.p, tl.n = tl.n,
                             p = p)

    scoresflat <- c(coord$score.mat)
    scoresnames <- paste0(rep(1:length(n), times = length(n)),
                          "-",
                          rep(1:length(n), each = length(n)))
    scoresnames <- scoresnames[!is.na(scoresflat)]
    scoresflat  <- scoresflat[!is.na(scoresflat)]
    names(scoresflat) <- scoresnames


    l = coord$U.T
    attr(l, "names") = "U.T"
    t = p
    attr(t, "names") = "p"
    result = list(statistic = l,
                  parameter = t,
                  coordinates = scoresflat,
                  method = "Data Driven k-Sample Umbrella Test")
    class(result) = c("htest", "ddst.test", "ddst.umbrella.test")

    result
  }
