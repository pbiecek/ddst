#' Data Driven Smooth Test for Normality; Unbounded Basis Functions
#'
#' Performs data driven smooth test for composite hypothesis of normality.
#' Null density is given by
#' \eqn{
#'   f(z;\gamma)=1/(\sqrt{2 \pi}\gamma_2) \exp(-(z-\gamma_1)^2/(2 \gamma_2^2))} for \eqn{z \in R}.
#'
#' @param x a (non-empty) numeric vector of data values
#' @param d.n an integer specifying the maximum dimension considered, only for advanced users
#' @param e.0 a (non-empty) numeric vector being the Monte Carlo estimate of the mean of the vector (C2; ... ;Cd.n) calculated using the function \code{ddst.normunbounded.bias()$e.0}
#' @param v.0 a (non-empty) numeric vector being the Monte Carlo estimate of the variance of the vector (C2; ... ;Cd.n) calculated using the function \code{ddst.normunbounded.bias()$e.0}
#' @param alpha a significance level
#' @param r.alpha a critical value of the alpha level R.n test
#' @param s.n.alpha a penalty in the auxiliary model selection rule
#' @param nr an integer specifying the number of runs for a critical value computation
#' @param compute.cv a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not
#' @param n sample size
#' @param ... further arguments
#'
#' @export
#'
#' @references
#' Ledwina, T., Wyłupek, G. (2015) Detection of non-Gaussianity by Ledwina and Wyłupek \emph{Journal of Statistical Computation and Simulation} 17, 3480-3497.
#'
#' @importFrom orthopolynom hermite.h.polynomials
#' @keywords htest
#' @examples
#' set.seed(7)
#' # H0 is true
#' z <- rnorm(100)
#' # let's look on first 20 coordinates
#' d.n <- 20
#'
#' # calculate finite sample corrections
#' # see 6.2. Composite null hypothesis H in the appendix materials
#' e.v <- ddst.normunbounded.bias(n = length(z))
#' e.v
#'
#' # simulated 1-alpha qunatiles, s(n, alpha) and  s.o(n, alpha)
#' # see Table 1 in the JSCS article
#' s <- 4.4
#' r.alpha <- 2.708
#'
#' t <- ddst.normubounded.test(z, d.n, e.v$e.0, e.v$v.0, r.alpha, s)
#' t
#' plot(t)
#'
#' # for Tephra data
#' z <- c(-1.748789, -1.75753, -1.740102, -1.740102, -1.731467, -1.765523,
#'        -1.761521, -1.72522, -1.80371, -1.745624, -1.872957, -1.729121,
#'        -1.81529, -1.888637, -1.887761, -1.881645, -1.91518, -1.849769,
#'        -1.755141, -1.665687, -1.764721, -1.736171, -1.736956, -1.737742,
#'        -1.687537, -1.804534, -1.790593, -1.808661, -1.784081, -1.729903,
#'        -1.711263, -1.748789, -1.772755, -1.72756, -1.71358, -1.821116,
#'        -1.839588, -1.839588, -1.830321, -1.807835, -1.747206, -1.788147,
#'        -1.759923, -1.786519, -1.726779, -1.738528, -1.754345, -1.781646,
#'        -1.641949, -1.755936, -1.775175, -1.736956, -1.705103, -1.743255,
#'        -1.82613, -1.826967, -1.780025, -1.684504, -1.751168)
#'
#' # calculate finite sample corrections
#' e.v <- ddst.normunbounded.bias(n = length(z))
#' e.v
#'
#' s <- 3.3
#' so <- 3.6
#' r.alpha <- 2.142
#'
#' t <- ddst.normubounded.test(z, d.n, e.v$e.0, e.v$v.0, r.alpha, s)
#' t
#' plot(t)
#' @rdname ddst.normunbounded.test
#'
`ddst.normubounded.test` <-
function(x,
         d.n = 20,
         e.0, v.0,
         r.alpha, s.n.alpha, alpha = 0.05,
         nr = 10000, compute.cv = TRUE) {
  n <- length(x)

  cal.C.unbounded.vec = (C.unbounded.vec(x, n, d.n) - e.0) / v.0
  cal.C.unbounded.vec.sq = cal.C.unbounded.vec ^ 2
  cal.N.d.n = cumsum(cal.C.unbounded.vec.sq)
  d.n.o = ceiling(d.n / 2)
  cal.Co.vec.sq = cal.C.unbounded.vec.sq[2 * 1:d.n.o]
  cal.No.d.n = cumsum(cal.Co.vec.sq)
  dim = 1:d.n
  dim.o = 1:d.n.o
  cal.A.s = which.max(cal.N.d.n - dim * s.n.alpha)
  cal.Ao.so = which.max(cal.No.d.n - dim.o * s.n.alpha)
  cal.A.2 = which.max(cal.N.d.n - dim * 2)
  cal.Ao.2 = which.max(cal.No.d.n - dim.o * 2)
  r.n = cal.unbounded.R(x, n)
  if (r.n <= r.alpha) {
    cal.A = cal.A.s
    cal.Ao = cal.Ao.so
  } else {
    cal.A = cal.A.2
    cal.Ao = cal.Ao.2
  }
  cal.N.cal.A.s = cal.N.d.n[cal.A.s]
  cal.N.cal.A = cal.N.d.n[cal.A]
  cal.N.cal.A.2 = cal.N.d.n[cal.A.2]
  cal.No.d.n.cal.Ao.so = cal.No.d.n[cal.Ao.so]
  cal.No.d.n.cal.Ao = cal.No.d.n[cal.Ao]
  cal.No.d.n.cal.Ao.2 = cal.No.d.n[cal.Ao.2]
  result = list(
    N.cal.A.s = cal.N.cal.A.s,
    N.cal.A = cal.N.cal.A,
    N.cal.A.2 = cal.N.cal.A.2,
    A.s = cal.A.s,
    A = cal.A,
    A.2 = cal.A.2,
    No.d.n.cal.Ao.so = cal.No.d.n.cal.Ao.so,
    No.d.n.cal.Ao = cal.No.d.n.cal.Ao,
    No.d.n.cal.Ao.2 = cal.No.d.n.cal.Ao.2,
    Ao.so = cal.Ao.so,
    Ao = cal.Ao,
    Ao.2 = cal.Ao.2,
    C.unbounded.vec = cal.C.unbounded.vec
  )


  l = result$A.s
  attr(l, "names") = "As"
  t = result$N.cal.A.s
  attr(t, "names") = "NAs"
  result = list(statistic = t,
                parameter = l,
                coordinates = result$C.unbounded.vec,
                method = "Data Driven Smooth Test for Normality - Unbounded Basis Functions")
  result$data.name = paste(paste(as.character(substitute(x)), collapse = ""),
                           ", d.n: ",
                           d.n,
                           ", r.alpha: ",
                           r.alpha,
                           ", s.n.alpha: ",
                           s.n.alpha,
                           sep = "")
  class(result) = c("htest", "ddst.test", "ddst.unbounded.test")
  return(result)
}

#' @export
#' @rdname ddst.normunbounded.test
`ddst.normunbounded.bias` <-
function(n = 100,
         d.n = 20,
         nr = 10000) {
  C.0 = matrix(0, d.n, nr)
  e.0 = matrix(0, 1, d.n)
  v.0 = matrix(0, 1, d.n)

  for (i in 1:nr) {
    x = rnorm(n)
    C.0[, i] = C.unbounded.vec(x, n, d.n)
  }
  for (j in 1:d.n) {
    e.0[1, j] = mean(C.0[j, ])
    v.0[1, j] = sd(C.0[j, ]) * sqrt((n - 1) / n)
  }
  list(e.0 = e.0,
       v.0 = v.0)
}

cal.unbounded.R = function(x, n) {
  S.sq = (sum((x - mean(x)) ^ 2)) / n
  score = n * (1 - delta.unbounded.k(x, n, 0) ^ 2 / S.sq)
  return(score)
}

# prepare hermite standardized
H.j.stan <- function(x, j) {
  fpol = as.function(pols.unbounded.max[[j + 1]])
  2 ^ (-j / 2) * fpol(x / sqrt(2)) / sqrt(factorial(j))
}
# prepare delta.k
delta.unbounded.k = function(x, n, k) {
  x.s = sort(x)
  ps = 0:n / n
  phi = dnorm(qnorm(ps))
  h.j = H.j.stan(qnorm(ps), k)
  h.j[1] = 0
  h.j[n + 1] = 0
  w.k = matrix(0, 1, n)
  for (j in 1:n) {
    w.k[1, j] = phi[j] * h.j[j] - phi[j + 1] * h.j[j + 1]
  }
  score = (w.k[1,] %*% x.s) / sqrt(k + 1)
  return(score)
}
# calculate C statistics
C.unbounded.vec = function(x, n, d.n) {
  score = matrix(0, 1, d.n)
  for (j in 1:d.n)
    score[1, j] = delta.unbounded.k(x, n, j) * sqrt((j + 2) * n) / delta.unbounded.k(x, n, 0)
  return(score[1, ])
}

# prepare hermite polynomials
d.n.unbounded.max <- 100
pols.unbounded.max <- hermite.h.polynomials(d.n.unbounded.max)
