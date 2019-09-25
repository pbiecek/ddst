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
        } else{
          W_T[l] = sum(B[l, 1:A[l]] ^ 2)
        }
      }
      score = sum((1 - p) * W_T)
      return(list(score = score, W_T = W_T, B = B))
    }

    test_W_T(x = x,
             n = n,
             d_N = d_N,
             c = c)
  }
