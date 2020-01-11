#' Data Driven Smooth Test Against Stochastic Dominance
#'
#' Performs data driven smooth non-parametric two-sample test against one-sided alternatives
#' (stochastic dominance).
#' Suppose that we have random samples from two distributions F and G.
#' The null hypothesis is that F(x) < G(x) for some x while the alternative is that at
#' F(x) >= G(x) for all x with strict inequality for at least one x.
#' Detailed description of the test statistic is provided in Ledwina and Wylupek (2012).
#'
#' @param x a (non-empty) numeric vector of data
#' @param y a (non-empty) numeric vector of data
#' @param k.n an integer specifying a level of complexity of the grid considered, only for advanced users
#' @param alpha a significance level
#' @param t an alpha-dependent tunning parameter in the penalty in the model selection rule
#' @param B an integer specifying the number of runs for a p-value and a critical value computation if any
#' @param compute.cv a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not
#'
#' @references Two-sample test against one-sided alternatives. Ledwina and Wylupek (2012). \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9469.2011.00787.x}
#' @export
#' @examples
#' set.seed(7)
#' # H0 is true
#' x <- runif(80)
#' y <- runif(80)
#' t <- ddst.againststochdom.test(x, y, alpha = 0.05, t = 2.2, k.n = 4)
#' t
#' plot(t)
#'
#' # H0 is false
#' # known fixed alternative
#' x <- runif(80)
#' y <- rbeta(80,4,2)
#' t <- ddst.againststochdom.test(x, y, alpha = 0.05, t = 2.2, k.n = 4)
#' t
#' plot(t)
#' @keywords htest
`ddst.againststochdom.test` <-
  function(x,
           y,
           k.n = 4,
           alpha = 0.05,
           t,# = 2.2,
           B = 10000,
           compute.cv = FALSE) {
    coord = ddst.stochasticdominance.Nk(x, y, t = t, k.N = k.n, alpha = alpha)    # coord square times n

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
                             ", k.n: ",
                             k.n,
                             sep = "")
    class(result) = c("htest", "ddst.test", "ddst.stochasticdominance.test")

    result
  }







# k-Sample Test
# Based on R code by Grzegorz Wylupek
# Two-sample test against one-sided alternatives.
# Ledwina and Wylupek (2012).
# https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9469.2011.00787.x
`ddst.stochasticdominance.Nk` <-
  function(x, y, t = 2.2, k.N = 4, alpha = 0.05) {
    ax.f <- function(k, kappa, varsigma, i) {
      2 ^ (-kappa - k - 2) * sqrt((2 * i - 1) * (2 * varsigma - 1) *
                                    (2 ^ (kappa + 1) - (2 * i - 1)) * (2 ^ (k +  1) - (2 * varsigma - 1)))
    }

    d.N = 2 ^ (k.N + 1) - 1

    Inv.Sigma = array(0, dim = c(d.N, d.N, k.N))
    for (k in 1:k.N) {
      for (kappa in 0:k) {
        for (i in 1:(2 ^ kappa)) {
          u <- 2 ^ kappa - 1 + i
          Inv.Sigma[u, u, k] = ax.f(kappa, kappa, i, i) * 2 ^ (k + 2)
        }
      }
      for (kappa in 0:(k - 1)) {
        for (i in 1:(2 ^ kappa)) {
          j <- 2 ^ kappa - 1 + i
          varsigma = 2 ^ (k - 1 - kappa) + (i - 1) * 2 ^ (k - kappa)
          p <- 2 ^ k - 1 + varsigma
          Inv.Sigma[j, p, k] <- -ax.f(k, kappa, varsigma, i) * 2 ^ (k + 1)
          Inv.Sigma[j, p + 1, k] <- -ax.f(k, kappa, varsigma + 1, i) * 2 ^ (k + 1)
        }

      }
      for (i in 1:(2 ^ k)) {
        u <- 2 ^ k - 1 + i
        for (s in 0:(k - 1)) {
          for (r in 1:(2 ^ s)) {
            v <- 2 ^ s - 1 + r
            Inv.Sigma[u, v, k] <- Inv.Sigma[v, u, k]
          }
        }
      }
    }
    l.j = function(k, i, z) {
      a.j = (2 * i - 1) / 2 ^ (k + 1)
      left = -sqrt((1 - a.j) / a.j)
      right = sqrt(a.j / (1 - a.j))
      score = ifelse(z < a.j, left, right)
      return(score)
    }
    L.j = function(k, i, x, y, m, n) {
      N = m + n
      H = ecdf(c(x, y))
      rx = H(x) - 1 / (2 * N)
      ry = H(y) - 1 / (2 * N)
      cx = sum(l.j(k, i, rx)) / m
      cy = sum(l.j(k, i, ry)) / n
      score = (cy - cx) * sqrt(m * n / N)
      return(score)
    }
    # x vector of the length m
    # y vector of the length n
    # H0: F=G | F<= G
    # HA: F >= G

    # Haar basis - 2i −1)/2k+1, k=0, 1, . . ., i=1, 2, . . ., 2k
    # t - some positive constant
    test.V = function(x, y, m, n, t, alpha) {
      N = m + n
      L.vec = matrix(0, 1, d.N)
      for (k in 0:k.N)
      {
        for (i in 1:(2 ^ k))
        {
          j = 2 ^ k - 1 + i
          L.vec[1, j] = L.j(k, i, x, y, m, n)
        }
      }
      W.d.vec = matrix(0, 1, k.N + 1)
      D.vec = 2 ^ (0:k.N + 1) - 1
      W.d.vec[1, 1] = L.vec[1, 1] ^ 2
      for (k in 2:(k.N + 1)) {
        W.d.vec[1, k] = L.vec[1, 1:D.vec[k]] %*% Inv.Sigma[1:D.vec[k], 1:D.vec[k], k - 1] %*% L.vec[1, 1:D.vec[k]]
      }
      S = which.max(W.d.vec[1, ] - D.vec * log(N))
      M = which.max(W.d.vec[1, ])

      if (isTRUE(max(abs(L.vec[1, ])) <= sqrt(t * log(N)))) {
        T = S
      } else{
        T = M
      }
      z.alpha = qnorm(alpha)
      V.S = sign(min(L.vec[1, ] - z.alpha)) * W.d.vec[1, S]
      V.T = sign(min(L.vec[1, ] - z.alpha)) * W.d.vec[1, T]

      return(list(V.T = V.T, V.S = V.S, T = T, S = S, L = L.vec[1, ]))
    }

    #  x <- runif(20)
    #  y <- runif(20)
    #  t <- 2.2
    m <- length(x)
    n <- length(y)

    test.V(x, y, m, n, t, alpha)

  }

