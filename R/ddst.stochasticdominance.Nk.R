# code based on Grzegorz Wylupek R Code
# Two-Sample Test Against One-Sided Alternatives
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

    result = c(V.T = V.T, V.S = V.S, T = T, S = S, L.vec[1, ])
    return(result)
  }

#  x <- runif(20)
#  y <- runif(20)
  #  t <- 2.2
  m <- length(x)
  n <- length(y)

  # QUESTION: TODO: Is V.T the test statistic?
  test.V(x, y, m, n, t, alpha)

}

