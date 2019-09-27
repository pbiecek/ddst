#' Data Driven Smooth Test for Uniformity
#'
#' Performs data driven smooth tests for simple hypothesis of uniformity on [0,1].
#' Embeding null model into the original exponential family introduced by Neyman (1937).
#' For more details see: \url{http://www.biecek.pl/R/ddst/description.pdf}.
#'
#' @aliases ddst.uniform.Nk
#' @param x a (non-empty) numeric vector of data values
#' @param base a function which returns orthogonal system, might be \code{ddst.base.legendre} for Legendre polynomials or \code{ddst.base.cos} for cosine system, see package description
#' @param c a parameter for model selection rule, see package description
#' @param B an integer specifying the number of replicates used in p-value computation
#' @param compute.p a logical value indicating whether to compute a p-value
#' @param Dmax an integer specifying the maximum number of coordinates, only for advanced users
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
#' @importFrom stats ecdf
#' @export
#'
#' @examples
#' # H0 is true
#' z <- runif(80)
#' t <- ddst.uniform.test(z, compute.p=TRUE)
#' t
#' plot(t)
#'
#' # known fixed alternative
#' z <- rnorm(80,10,16)
#' t <- ddst.uniform.test(pnorm(z, 10, 16), compute.p=TRUE)
#' t
#' plot(t)
#'
#' # H0 is false
#' z <- rbeta(80,4,2)
#' (t <- ddst.uniform.test(z, compute.p=TRUE))
#' t$p.value
#' plot(t)
#' @keywords htest
`ddst.uniform.test` <-
  function(x,
           base = ddst.base.legendre,
           c = 2.4,
           B = 1000,
           compute.p = FALSE,
           Dmax = 10,
           ...) {
    # method.name = as.character(substitute(base))
    # only Legendre is implemented yet
    base = ddst.base.legendre
    method.name = "ddst.base.legendre"

    n = length(x)
    if (n < 5)
      stop("length(x) should be at least 5")
    coord = ddst.uniform.Nk(x, base, Dmax = Dmax)    # coord square times n
    l = ddst.IIC(coord, n, c)
    attr(l, "names") = "n. coord"
    t = coord[l]
    attr(t, "names") = "WT"
    result = list(statistic = t,
                  parameter = l,
                  coordinates = coord - c(0, coord[-Dmax]),
                  method = "Data Driven Smooth Test for Uniformity")
    result$data.name = paste(paste(as.character(substitute(x)), collapse =
                                     ""),
                             ",   base: ",
                             method.name,
                             "   c: ",
                             c,
                             sep = "")
    class(result) = c("htest", "ddst.test")
    if (compute.p) {
      tmp = numeric(B)
      for (i in 1:B) {
        y = runif(n)
        tmpC = ddst.uniform.Nk(y, base, Dmax = Dmax)
        l = ddst.IIC(tmpC, n, c)
        tmp[i] = tmpC[l]
      }
      result$p.value = mean(tmp > t)
    }
    result
  }
