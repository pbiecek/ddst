#' Data Driven Smooth Test for Umbrella Alternatives; Known Peak
#'
#' Performs data driven smooth test for so-called umbrella alternatives in k-sample problem.
#' Suppose that we have random samples from k distributions F_i where i = 1, ..., k.
#' The null hypothesis is that there is no umbrella pattern,
#' i.e. F_1 >= ... >= F_p <= ... <= F_k and F_i != F_j for some i and j.
#' The alternative is that there is an umbrella pattern
#' i.e. F_1 >= ... >= F_p <= ... <= F_k and F_i != F_j for some i and j.
#' Detailed description of the test statistic is provided in Wylupek (2016).
#'
#' @param x a list of k (non-empty) numeric vectors of data
#' @param p a peak of the umbrella pattern
#' @param r.N a (p(p-1)=2+(k-p)(k-p+1)/2)-dimensional vector specifying the levels of complexity of the grids considered, only for advanced users
#' @param alpha a significance level
#' @param t.p an alpha-dependent (p(p-1)=2+(k-p)(k-p+1)/2)-dimensional vector of the tunning parameters in the penalties in the model selection rules T.o
#' @param t.n an alpha-dependent (k-1)-dimensional vector of the tunning parameters in the penalties in the model selection rules T.tilde
#' @param nr an integer specifying the number of runs for a p-value and a critical value computation if any
#' @param compute.cv a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not
#'
#' @references An automatic test for the umbrella alternatives. Wylupek (2016) \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/sjos.12231}
#' @export
#' @examples
#' set.seed(7)
#' # H0 is true
#' x = runif(80)
#' y = runif(80) + 0.2
#' z = runif(80)
#' t <- ddst.umbrellaknownp.test(list(x, y, z), p = 2, t.p = 2.2, t.n = 2.2)
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
#' t <- ddst.umbrellaknownp.test(list(x1, x2, x3, x4, x5, x6, x7), p = 3, t.p = 2.2, t.n = 2.2)
#' t
#' plot(t)
#'
#' t <- ddst.umbrellaknownp.test(list(x1, x2, x3, x4, x5, x6, x7), p = 5, t.p = 2.2, t.n = 2.2)
#' t
#' plot(t)
#'
#' @keywords htest
`ddst.umbrellaknownp.test` <-
  function(x,
           p,
           r.N  = rep(4, length(x) - 1),
           alpha = 0.05,
           t.p, t.n,
           nr = 100000, compute.cv = TRUE) {
##    if (is.list(x)) {
## now it has to be a list
      # x is list with coordinates
      x.vector <- unlist(x)
      n <- sapply(x, length)
##    } else {
##      # x is a vector, other vectors are in ...
##      x.vector <- unlist(c(x, list(...)))
##      n <- c(length(x),
##             sapply(list(...), length))
##    }

    coord = ddst.umbrella.Nk(x.vector, n,
                             tlh.p = t.p, tl.n = t.n,
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
    result$data.name = paste(paste(as.character(substitute(x)), collapse = ""),
                             ", t.p: ",
                             t.p,
                             ", t.n: ",
                             t.n,
                             ", p: ",
                             p,
                             sep = "")

    result
  }









# Umbrella Test
# Based on R code by Grzegorz Wylupek
# An automatic test for the umbrella alternatives.
# Wylupek (2016)
# https://onlinelibrary.wiley.com/doi/abs/10.1111/sjos.12231
# sjos12231-sup-0002-supplementary.r
`ddst.umbrella.Nk` <-
  function(x, n, tlh.p = 2.2, tl.n = 2.2, p = 3) {
    # what about default for tlh.p

    l.j = function(r, i, v) {
      a.j = (2 * i - 1) / 2 ^ (r + 1)
      left = -sqrt((1 - a.j) / a.j)
      right = sqrt(a.j / (1 - a.j))
      score = ifelse(v < a.j, left, right)
      return(score)
    }

    L.j.lh = function(r, i, xl, xh, nl, nh) {
      Nlh = nl + nh
      Hlh = ecdf(c(xl, xh))
      rxl = Hlh(xl) - 1 / (2 * Nlh)
      rxh = Hlh(xh) - 1 / (2 * Nlh)
      cxl = mean(l.j(r, i, rxl))
      cxh = mean(l.j(r, i, rxh))
      score = (cxh - cxl) * sqrt(nl * nh / Nlh)
      return(score)
    }

    r.Nlh = 4
    d.Nlh = 2 ^ (r.Nlh + 1) - 1

    Z.l = function(xl, xl1, nl, nl1, tl.n, alpha.l) {
      Nl = nl + nl1
      L.vec.l = matrix(0, 1, d.Nlh)
      L.vec.l.n = matrix(0, 1, d.Nlh)
      for (r in 0:r.Nlh) {
        for (i in 1:(2 ^ r)) {
          j = 2 ^ r - 1 + i
          L.vec.l[1, j] = L.j.lh(r, i, xl, xl1, nl, nl1)
          L.vec.l.n[1, j] = max(-L.vec.l[1, j], 0)
        }
      }
      Q.d.vec.l.n = matrix(0, 1, r.Nlh + 1)
      D.vec.l = 2 ^ (0:r.Nlh + 1) - 1
      for (r in 1:(r.Nlh + 1)) {
        Q.d.vec.l.n[1, r] = L.vec.l.n[1, 1:D.vec.l[r]] %*% L.vec.l.n[1, 1:D.vec.l[r]]
      }
      Sl.n = which.max(Q.d.vec.l.n[1, ] - D.vec.l * log(Nl))
      Ml.n = which.max(Q.d.vec.l.n[1, ])
      if (max(L.vec.l.n[1, ]) <= sqrt(tl.n * log(Nl))) {
        Tl.n = Sl.n
      } else{
        Tl.n = Ml.n
      }
      z.alpha.l = qnorm(alpha.l)
      Z.Sl.n = sign(min(L.vec.l[1, 1:D.vec.l[Sl.n]] - z.alpha.l))
      Z.Tl.n = sign(min(L.vec.l[1, 1:D.vec.l[Tl.n]] - z.alpha.l))
      Z.Ml.n = sign(min(L.vec.l[1, 1:D.vec.l[Ml.n]] - z.alpha.l))
      result = c(Z.Sl.n, Z.Tl.n, Z.Ml.n, Sl.n, Tl.n, Ml.n)
      return(result)
    }

    Q.lh = function(xl, xh, nl, nh, tlh.p) {
      Nlh = nl + nh
      L.vec.lh = matrix(0, 1, d.Nlh)
      L.vec.lh.p = matrix(0, 1, d.Nlh)
      for (r in 0:r.Nlh) {
        for (i in 1:(2 ^ r)) {
          j = 2 ^ r - 1 + i
          L.vec.lh[1, j] = L.j.lh(r, i, xl, xh, nl, nh)
          L.vec.lh.p[1, j] = max(L.vec.lh[1, j], 0)
        }
      }
      Q.d.vec.lh.p = matrix(0, 1, r.Nlh + 1)
      D.vec.lh = 2 ^ (0:r.Nlh + 1) - 1
      for (r in 1:(r.Nlh + 1)) {
        Q.d.vec.lh.p[1, r] = L.vec.lh.p[1, 1:D.vec.lh[r]] %*% L.vec.lh.p[1, 1:D.vec.lh[r]]
      }
      Slh.p = which.max(Q.d.vec.lh.p[1, ] - D.vec.lh * log(Nlh))
      Mlh.p = which.max(Q.d.vec.lh.p[1, ])
      if (max(L.vec.lh.p[1, ]) <= sqrt(tlh.p * log(Nlh))) {
        Tlh.p = Slh.p
      } else{
        Tlh.p = Mlh.p
      }
      Q.Slh.p = Q.d.vec.lh.p[1, Slh.p]
      Q.Tlh.p = Q.d.vec.lh.p[1, Tlh.p]
      Q.Mlh.p = Q.d.vec.lh.p[1, Mlh.p]
      result = c(Q.Slh.p, Q.Tlh.p, Q.Mlh.p, Slh.p, Tlh.p, Mlh.p)
      return(result)
    }

    test.U = function(x, n, tlh.p, tl.n, alpha = 0.05, p = 3) {
      k = length(n)
      N = sum(n)
      nc = c(0, cumsum(n))
      p.mat = matrix(0, k, k)
      score.mat = matrix(NA, k, k)
      for (l in 1:k) {
        for (h in 1:k) {
          p.mat[l, h] = (n[l] + n[h]) / N
        }
      }
      p.s = 0
      for (j in 1:(k - 1)) {
        p.s = p.s + p.mat[j, j + 1]
      }
      w = c()
      for (l in 1:(k - 1)) {
        w[l] = p.mat[l, l + 1] / p.s
      }
      alpha.vec = c()
      alpha.vec = w * alpha
      Q.S.o = 0
      Q.T.o = 0
      Q.M.o = 0
      S.p.vec = c()
      T.p.vec = c()
      M.p.vec = c()
      Z.S.t.vec = c()
      Z.T.t.vec = c()
      Z.M.t.vec = c()
      S.n.vec = c()
      T.n.vec = c()
      M.n.vec = c()
      ind = 0
      if (p == 1) {
        for (l in 1:(k - 1)) {
          for (h in (l + 1):k) {
            score = Q.lh(x[(nc[h] + 1):nc[h + 1]], x[(nc[l] + 1):nc[l + 1]], n[h], n[l], tlh.p)
            score.mat[l,h] <- score[2]
            Q.S.o = Q.S.o + p.mat[h, l] * score[1]
            Q.T.o = Q.T.o + p.mat[h, l] * score[2]
            Q.M.o = Q.M.o + p.mat[h, l] * score[3]
            ind = ind + 1
            S.p.vec[ind] = score[4]
            T.p.vec[ind] = score[5]
            M.p.vec[ind] = score[6]
          }
        }
        for (l in 1:(k - 1)) {
          score = Z.l(x[(nc[l + 1] + 1):nc[l + 2]], x[(nc[l] + 1):nc[l + 1]], n[l +
                                                                                  1], n[l], tl.n, alpha.vec[l])
          Z.S.t.vec[l] = score[1]
          Z.T.t.vec[l] = score[2]
          Z.M.t.vec[l] = score[3]
          S.n.vec[l] = score[4]
          T.n.vec[l] = score[5]
          M.n.vec[l] = score[6]
        }
        U.S = min(Z.S.t.vec) * Q.S.o
        U.T = min(Z.T.t.vec) * Q.T.o
        U.M = min(Z.M.t.vec) * Q.M.o
        result = c(U.S,
                   U.T,
                   U.M,
                   S.p.vec,
                   T.p.vec,
                   M.p.vec,
                   S.n.vec,
                   T.n.vec,
                   M.n.vec)
      }
      if (p == k) {
        for (l in 1:(k - 1)) {
          for (h in (l + 1):k) {
            score = Q.lh(x[(nc[l] + 1):nc[l + 1]], x[(nc[h] + 1):nc[h + 1]], n[l], n[h], tlh.p)
            score.mat[l,h] <- score[2]
            Q.S.o = Q.S.o + p.mat[l, h] * score[1]
            Q.T.o = Q.T.o + p.mat[l, h] * score[2]
            Q.M.o = Q.M.o + p.mat[l, h] * score[3]
            ind = ind + 1
            S.p.vec[ind] = score[4]
            T.p.vec[ind] = score[5]
            M.p.vec[ind] = score[6]
          }
        }
        for (l in 1:(k - 1)) {
          score = Z.l(x[(nc[l] + 1):nc[l + 1]], x[(nc[l + 1] + 1):nc[l + 2]], n[l], n[l +
                                                                                        1], tl.n, alpha.vec[l])
          Z.S.t.vec[l] = score[1]
          Z.T.t.vec[l] = score[2]
          Z.M.t.vec[l] = score[3]
          S.n.vec[l] = score[4]
          T.n.vec[l] = score[5]
          M.n.vec[l] = score[6]
        }
        U.S = min(Z.S.t.vec) * Q.S.o
        U.T = min(Z.T.t.vec) * Q.T.o
        U.M = min(Z.M.t.vec) * Q.M.o
        result = c(U.S,
                   U.T,
                   U.M,
                   S.p.vec,
                   T.p.vec,
                   M.p.vec,
                   S.n.vec,
                   T.n.vec,
                   M.n.vec)
      }
      if (p != 1 & p != k) {
        for (l in 1:(p - 1)) {
          for (h in (l + 1):p) {
            score = Q.lh(x[(nc[l] + 1):nc[l + 1]], x[(nc[h] + 1):nc[h + 1]], n[l], n[h], tlh.p)
            score.mat[l,h] <- score[2]
            Q.S.o = Q.S.o + p.mat[l, h] * score[1]
            Q.T.o = Q.T.o + p.mat[l, h] * score[2]
            Q.M.o = Q.M.o + p.mat[l, h] * score[3]
            ind = ind + 1
            S.p.vec[ind] = score[4]
            T.p.vec[ind] = score[5]
            M.p.vec[ind] = score[6]
          }
        }
        for (l in p:(k - 1)) {
          for (h in (l + 1):k) {
            score = Q.lh(x[(nc[h] + 1):nc[h + 1]], x[(nc[l] + 1):nc[l + 1]], n[h], n[l], tlh.p)
            score.mat[l,h] <- score[2]
            Q.S.o = Q.S.o + p.mat[h, l] * score[1]
            Q.T.o = Q.T.o + p.mat[h, l] * score[2]
            Q.M.o = Q.M.o + p.mat[h, l] * score[3]
            ind = ind + 1
            S.p.vec[ind] = score[4]
            T.p.vec[ind] = score[5]
            M.p.vec[ind] = score[6]
          }
        }
        for (l in 1:(p - 1)) {
          score = Z.l(x[(nc[l] + 1):nc[l + 1]], x[(nc[l + 1] + 1):nc[l + 2]], n[l], n[l +
                                                                                        1], tl.n, alpha.vec[l])
          Z.S.t.vec[l] = score[1]
          Z.T.t.vec[l] = score[2]
          Z.M.t.vec[l] = score[3]
          S.n.vec[l] = score[4]
          T.n.vec[l] = score[5]
          M.n.vec[l] = score[6]
        }
        for (l in p:(k - 1)) {
          score = Z.l(x[(nc[l + 1] + 1):nc[l + 2]], x[(nc[l] + 1):nc[l + 1]], n[l +
                                                                                  1], n[l], tl.n, alpha.vec[l])
          Z.S.t.vec[l] = score[1]
          Z.T.t.vec[l] = score[2]
          Z.M.t.vec[l] = score[3]
          S.n.vec[l] = score[4]
          T.n.vec[l] = score[5]
          M.n.vec[l] = score[6]
        }
        U.S = min(Z.S.t.vec) * Q.S.o
        U.T = min(Z.T.t.vec) * Q.T.o
        U.M = min(Z.M.t.vec) * Q.M.o
        result = c(U.S,
                   U.T,
                   U.M,
                   S.p.vec,
                   T.p.vec,
                   M.p.vec,
                   S.n.vec,
                   T.n.vec,
                   M.n.vec)
      }
      return(list(U.T = U.T, T.p.vec = T.p.vec, T.n.vec = T.n.vec, score.mat = score.mat))
    }

    # select UT - US - UM
    test.U(x, n, tlh.p = tlh.p, tl.n = tl.n, p = p)
  }
