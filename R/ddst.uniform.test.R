#' Data Driven Smooth Test for Uniformity
#'
#' Performs data driven smooth tests for simple hypothesis of uniformity on [0,1].
#' Embeding null model into the original exponential family introduced by Neyman (1937).
#'
#' @aliases ddst.uniform.Nk
#' @param x a (non-empty) numeric vector of data
#' @param base a function which returns an orthonormal system, possible choice: \code{ddst.base.legendre} for the Legendre polynomials and \code{ddst.base.cos} for the cosine system
#' @param d.n an integer specifying the maximum dimension considered, only for advanced users
#' @param c a calibrating parameter in the penalty in the model selection rule
#' @param nr an integer specifying the number of runs for a p-value and a critical value computation if any
#' @param compute.p a logical value indicating whether to compute a p-value or not
#' @param alpha a significance level
#' @param compute.cv a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not
#' @param ... further arguments
#'
#' @return
#' An object of class \code{htest}
#' \item{statistic }{the value of the test statistic.}
#' \item{parameter }{the number of choosen coordinates (k).}
#' \item{method }{a character string indicating the parameters of performed test. }
#' \item{data.name }{a character string giving the name(s) of the data. }
#' \item{p.value }{the p-value for the test, computed only if \code{compute.p=T}.}
#'
#' @references
#' Inglot, T., Ledwina, T. (2006). Towards data driven selection of a penalty function for data driven Neyman tests. \emph{ Linear Algebra and its Appl.} \bold{ 417}, 579--590.
#'
#' Ledwina, T. (1994). Data driven version of Neyman's smooth test of fit. \emph{ J. Amer. Statist. Assoc.} \bold{ 89} 1000-1005.
#'
#' Neyman, J. (1937). `Smooth test' for goodness of fit. \emph{Skand. Aktuarietidskr.} \bold{ 20}, 149-199.
#'
#' @importFrom stats ecdf dnorm quantile sd
#' @export
#'
#' @examples
#' set.seed(7)
#' # H0 is true
#' z <- runif(80)
#' \dontrun{
#' t <- ddst.uniform.test(z, compute.p = TRUE, d.n = 10)
#' t
#' plot(t)
#'
#' # known fixed alternative
#' z <- rnorm(80,10,16)
#' t <- ddst.uniform.test(pnorm(z, 10, 16), compute.p = TRUE, d.n = 10)
#' t
#' plot(t)
#'
#' # H0 is false
#' z <- rbeta(80,4,2)
#' (t <- ddst.uniform.test(z, compute.p = TRUE, d.n = 10))
#' t$p.value
#' plot(t)
#' }
#' @keywords htest
`ddst.uniform.test` <-
  function(x,
           base = ddst.base.legendre,
           d.n = 10,
           c = 2.4,
           nr = 100000,
           compute.p = TRUE,
           alpha = 0.05,
           compute.cv = TRUE,
           ...) {
    # method.name = as.character(substitute(base))
    # only Legendre is implemented yet
    base = ddst.base.legendre
    method.name = "ddst.base.legendre"

    n = length(x)
    if (n < 5)
      stop("length(x) should be at least 5")
    coord = ddst.uniform.Nk(x, base, Dmax = d.n)    # coord square times n
    l = ddst.IIC(coord, n, c)
    attr(l, "names") = "T"
    t = coord[l]
    attr(t, "names") = "WT"
    result = list(statistic = t,
                  parameter = l,
                  coordinates = coord - c(0, coord[-d.n]),
                  method = "Data Driven Smooth Test for Uniformity")

        if (compute.p | compute.cv) {
      tmp = numeric(nr)
      for (i in 1:nr) {
        y = runif(n)
        tmpC = ddst.uniform.Nk(y, base, Dmax = d.n)
        l = ddst.IIC(tmpC, n, c)
        tmp[i] = tmpC[l]
      }
      if (compute.p) {
        result$p.value = mean(tmp > t)
      }
      if (compute.cv) {
        result$cv = quantile(tmp, alpha)
      }
    }

    result$data.name = paste(paste(as.character(substitute(x)), collapse =
                                     ""),
                             ", base: ",
                             method.name,
                             "  c: ",
                             c,
                             "  d.n: ",
                             d.n,
                             ifelse(compute.cv, paste0("  cv(",alpha,") : ",signif(result$cv, 5),")"), ""),
                             sep = "")
    class(result) = c("htest", "ddst.test")

    result
  }






`ddst.uniform.Nk` <-
  function(x, base = ddst.base.legendre, Dmax = 10) {
    n = length(x)
    maxN = max(min(Dmax, n-2, 20),1)
    coord = numeric(maxN)
    for (j in 1:maxN)
      coord[j] = ddst.phi(x, j, base)
    coord = cumsum(coord^2*n)
  }


