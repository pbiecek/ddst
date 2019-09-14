#' Data Driven Smooth Test for Uniformity
#'
#' Performs data driven smooth tests for simple hypothesis of uniformity on [0,1].
#'
#' Embeding null model into the original exponential family introduced by Neyman (1937) leads to the information matrix \emph{ I} being identity and smooth test statistic with \emph{k} components
#' $$
#' W_k=1/\\sqrt(n) \\sum_\{j=1\}^k sum_\{i=1\}^n [\\phi_j(Z_i)]^2,
#' $$
#' where $phi_j$ is $j$th degree normalized Legendre polynomial on [0,1] (default value of parameter base = `'ddst.base.legendre'`). Alternatively, in our implementation, cosine system can be selected (base = `'ddst.base.cos'`). For details see Ledwina (1994) and Inglot and Ledwina (2006).
#'
#' An application of the pertaining selection rule \emph{T} for choosing \emph{k} gives related `ddst.uniform.test()' based on statistic \emph{$W_T$}.
#'
#' Similar approach applies to testing goodness-of-fit to any fully specified continuous distribution function \emph{F}. For this purpose it is enough to apply the above solution to transformed observations \emph{$F(z_1),...,F(z_n)$}.
#'
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
#' z = runif(80)
#' ddst.uniform.test(z, compute.p=TRUE)
#'
#' # known fixed alternative
#' z = rnorm(80,10,16)
#' ddst.uniform.test(pnorm(z, 10, 16), compute.p=TRUE)
#'
#' # H0 is false
#' z = rbeta(80,4,2)
#' (t = ddst.uniform.test(z, compute.p=TRUE))
#' t$p.value
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
                  method = "Data Driven Smooth Test for Uniformity")
    result$data.name = paste(paste(as.character(substitute(x)), collapse =
                                     ""),
                             ",   base: ",
                             method.name,
                             "   c: ",
                             c,
                             sep = "")
    class(result) = "htest"
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
