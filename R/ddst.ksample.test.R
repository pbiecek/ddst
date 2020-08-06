#' Data Driven Smooth Test for k-Sample Problem
#'
#' Performs data driven smooth test for the classical k-sample problem.
#' Suppose that we have random samples from k distributions F_i where i = 1, ..., k.
#' The null hypothesis is that F_1 = ... = F_k while the alternative is that at
#' least two distributions are different.
#' Detailed description of the test statistic is provided in Wylupek (2010).
#'
#' @param x a list of k (non-empty) numeric vectors of data
#' @param d.N an integer specifying the maximum dimension considered, only for advanced users
#' @param c a calibrating parameter in the penalty in the model selection rule
#' @param nr an integer specifying the number of runs for a p-value and a critical value computation if any
#' @param compute.p a logical value indicating whether to compute a p-value or not
#' @param alpha a significance level
#' @param compute.cv a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not
#'
#' @references Data-driven k-sample tests. Wylupek (2010) \url{https://www.jstor.org/stable/40586684?seq=1}
#' @export
#' @examples
#' set.seed(7)
#' # H0 is false
#' x <- runif(80)
#' y <- rexp(80, 1)
#' z <- runif(80)
#' t <- ddst.ksample.test(list(x, y, z))
#' t
#' plot(t)
#'
#' # H0 is true
#' x <- runif(80)
#' y <- runif(80)
#' z <- runif(80)
#' t <- ddst.ksample.test(list(x, y, z))
#' t
#' plot(t)
#'
#' @keywords htest
`ddst.ksample.test` <-
  function(x,
           d.N = 12,
           c = 2.3,
           nr = 100000,
           compute.p = TRUE,
           alpha = 0.05,
           compute.cv = TRUE) {
##    if (is.list(x)) {
      # x is list with coordinates
      x.vector <- unlist(x)
      n <- sapply(x, length)
##    } else {
##      # x is a vector, other vectors are in ...
##      x.vector <- unlist(c(x, list(...)))
##      n <- c(length(x),
##             sapply(list(...), length))
##    }
    coord = ddst.ksample.Nk(x.vector, n, d_N = d.n, c = c)

    l = coord$score
    attr(l, "names") = "WT"
    t = paste0("(", paste0(coord$L_T, collapse = ", "), ")")
    attr(t, "names") = "T"
    result = list(statistic = l,
                  parameter = t,
                  W.T.l = coord$W_T,
                  coordinates = coord$B,
                  method = "Data Driven Smooth Test for Two-Sample Problem")
    result$data.name = paste(paste(as.character(substitute(x)), collapse = ""),
                             ", c: ",
                             c,
                             ", d.n: ",
                             d.n,
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









# k-Sample Test
# Based on R code by Grzegorz Wylupek
# Data-driven k-sample tests.
# Wylupek (2010)
# https://www.jstor.org/stable/40586684?seq=1
`ddst.ksample.Nk` <-
  function(x, n, d_N = 12, c = 2.3) {
    pols <- slegendre.polynomials(d_N, normalized = TRUE)
    b_j <- function(u, j) {
      fp = as.function(pols[[j + 1]])
      fp(u)
    }
    Bl <- function(xl, xlc, nl, nlc, d_N) {
      N = nl + nlc
      H = ecdf(c(xl, xlc))
      rxl = H(xl) - 1 / (2 * N)
      rxlc = H(xlc) - 1 / (2 * N)
      cxl = NULL
      cxlc = NULL
      for (j in 1:d_N) {
        cxl[j] = mean(b_j(rxl, j))
        cxlc[j] = mean(b_j(rxlc, j))
      }
      sqrt(nl * nlc / N) * (cxl - cxlc)
    }
    test_W_T <- function(x, n, d_N = 12, c = 2.3) {
      N = sum(n)
      p = n / N
      k = length(n)
      n_cum = c(0, cumsum(n))
      B = matrix(0, k, d_N)
      S = NULL
      A = NULL
      W_T = NULL
      L_T = NULL
      for (l in 1:k) {
        nl = n[l]
        nlc = N - n[l]
        xl = x[(n_cum[l] + 1):n_cum[l + 1]]
        xlc = x[-((n_cum[l] + 1):n_cum[l + 1])]
        B[l, ] = Bl(xl, xlc, nl, nlc, d_N)
        S[l] = which.max(cumsum(B[l, ] ^ 2) - 1:d_N * log(N))
        A[l] = which.max(cumsum(B[l, ] ^ 2) - 1:d_N * 2)
        if (max(abs(B[l, ])) <= sqrt(c * log(N))) {
          W_T[l] = sum(B[l, 1:S[l]] ^ 2)
          L_T[l] = S[l]
        } else{
          W_T[l] = sum(B[l, 1:A[l]] ^ 2)
          L_T[l] = A[l]
        }
      }
      score = sum((1 - p) * W_T)
      return(list(score = score, W_T = W_T, B = B, L_T = L_T))
    }

    test_W_T(x = x,
             n = n,
             d_N = d_N,
             c = c)
  }

