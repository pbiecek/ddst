#' Data Driven Smooth Test for Umbrella Alternatives; Unknown Peak
#'
#' @param x a list of k (non-empty) numeric vectors of data
#' @param r.N a (p(p-1)=2+(k-p)(k-p+1)/2)-dimensional vector specifying the levels of complexity of the grids considered, only for advanced users
#' @param alpha a significance level
#' @param t.p.aux an auxiliary alpha-dependent k x (k - 1) matrix of the tunning parameters in the penalties in the model selection rules To aux employed for estimation of the peak
#' @param t.n an alpha-dependent (k-1)-dimensional vector of the tunning parameters in the penalties in the model selection rules T.tilde
#' @param t.p an alpha-dependent (p(p-1)=2+(k-p)(k-p+1)/2)-dimensional vector of the tunning parameters in the penalties in the model selection rules T.o
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
`ddst.umbrellaunknownp.test` <-
  function(x,
           r.N  = rep(4, length(x) - 1),
           alpha = 0.05,
           t.p.aux,
           t.p, t.n,
           nr = 100000, compute.cv = TRUE) {

    statistics <- sapply(seq_along(x), function(p)
      ddst.umbrellaknownp.test(x, p = p, r.N = r.N, t.p = t.p, t.n = t.n, nr = nr, compute.cv = compute.cv))

    max(statistics)
  }

