#' Data Driven k-Sample Test
#'
#' Performs data driven smooth test for the classical k-sample problem.
#' Suppose that we have random samples from k distributions F_i where i = 1, ..., k.
#' The null hypothesis is that F_1 = ... = F_k while the alternative is that at
#' least two distributions are different.
#' Detailed description of the test statistic is provided in Wylupek (2010).
#'
#' @param x a list with k-vectors or a single vector.
#' @param ... if x is a single vector, then remaing k-1 vectors are provided in the ... argument
#' @param d_N an integer, number of coordinates that measure potential deviation from null hypothesis
#' @param c a positive number, penalty for model selection rule. Section 5.1 in Wylupek (2010), suggest that a good choice is c = 2, when k = 2, and  c = 2.3, when k >= 3.
#'
#' @references Data-driven k-sample tests. Wylupek (2010) \url{https://www.jstor.org/stable/40586684?seq=1}
#' @export
#' @examples
#' # H0 is false
#' x <- runif(80)
#' y <- rexp(80, 1)
#' z <- runif(80)
#' t <- ddst.ksample.test(x, y, z)
#' t
#' plot(t)
#'
#' # H0 is true
#' x <- runif(80)
#' y <- runif(80)
#' z <- runif(80)
#' t <- ddst.ksample.test(x, y, z)
#' t <- ddst.ksample.test(list(x, y, z))
#' t
#' plot(t)
#'
#' @keywords htest
`ddst.ksample.test` <-
  function(x,
           ...,
           d_N = 12, c = 2.3) {
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

    l = coord$score
    attr(l, "names") = "WT"
    t = paste0("(", paste0(coord$L_T, collapse = ", "), ")")
    attr(t, "names") = "T"
    result = list(statistic = l,
                  parameter = t,
                  W.T.l = coord$W_T,
                  coordinates = coord$B,
                  method = "Data Driven k-Sample Test")
    result$data.name = paste(paste(as.character(substitute(x)), collapse = ""),
                             ", Dmax: ",
                             d_N,
                             sep = "")
    class(result) = c("htest", "ddst.test", "ddst.ksample.test")

    result
  }

# remove some wired error
polynomial <- function (coef = c(0, 1)) {
  a <- as.numeric(coef)
  while ((la <- length(a)) > 1 && a[la] == 0) a <- a[-la]
  structure(a, class = "polynomial")
}
