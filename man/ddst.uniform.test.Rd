\name{ddst.uniform.test}
\alias{ddst.uniform.test}
\alias{ddst.uniform.Nk}
\title{ Data Driven Smooth Test for Uniformity }
\description{
Performs data driven smooth tests for simple hypothesis of uniformity on [0,1].
}
\usage{
ddst.uniform.test(x, base = ddst.base.legendre, c = 2.4, B = 1000, compute.p = F,
    Dmax = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  a (non-empty) numeric vector of data values. }
  \item{base}{ a function which returns orthogonal system, might be \code{ddst.base.legendre} for Legendre polynomials or \code{ddst.base.cos} for cosine system, see package description. }
  \item{c}{ a parameter for model selection rule, see package description. }
  \item{B}{ an integer specifying the number of replicates used in p-value computation. }
  \item{compute.p}{  a logical value indicating whether to compute a p-value. }
  \item{Dmax}{ an integer specifying the maximum number of coordinates, only for advanced users. }
  \item{\dots}{ further arguments. }
}
\details{
Embeding null model into the original exponential family introduced by Neyman (1937) leads to the information matrix \emph{ I} being identity and smooth test statistic with \emph{k} components
\emph{
$W_k=[1/sqrt(n) sum_{j=1}^k sum_{i=1}^n phi_j(Z_i)]^2$},
where \emph{$phi_j$} is \emph{j}th degree normalized Legendre polynomial on [0,1] (default value of parameter base = `ddst.base.legendre'). Alternatively, in our implementation, cosine system can be selected (base = `ddst.base.cos'). For details see Ledwina (1994) and Inglot and Ledwina (2006).

An application of the pertaining selection rule \emph{T} for choosing \emph{k} gives related `ddst.uniform.test()' based on statistic \emph{$W_T$}.

Similar approach applies to testing goodness-of-fit to any fully specified continuous distribution function \emph{F}. For this purpose it is enough to apply the above solution to transformed observations \emph{$F(z_1),...,F(z_n)$}.

For more details see: \url{http://www.biecek.pl/R/ddst/description.pdf}.

}
\value{
  An object of class \code{htest}
  \item{statistic }{the value of the test statistic.}
  \item{parameter }{the number of choosen coordinates (k).}
  \item{method }{a character string indicating the parameters of performed test. }
  \item{data.name }{a character string giving the name(s) of the data. }
  \item{p.value }{the p-value for the test, computed only if \code{compute.p=T}.}
}
\references{
Inglot, T., Ledwina, T. (2006). Towards data driven selection of a penalty function for data driven Neyman tests. \emph{ Linear Algebra and its Appl.} \bold{ 417}, 579--590.

Ledwina, T. (1994). Data driven version of Neyman's smooth test of fit. \emph{ J. Amer. Statist. Assoc.} \bold{ 89} 1000-1005.

Neyman, J. (1937). `Smooth test' for goodness of fit. \emph{Skand. Aktuarietidskr.} \bold{ 20}, 149-199.
}
\author{ Przemyslaw Biecek and Teresa Ledwina }
\examples{

# H0 is true
z = runif(80)
ddst.uniform.test(z, compute.p=TRUE)

# known fixed alternative
z = rnorm(80,10,16)
ddst.uniform.test(pnorm(z, 10, 16), compute.p=TRUE)


# H0 is false
z = rbeta(80,4,2)
(t = ddst.uniform.test(z, compute.p=TRUE))
t$p.value

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{nonparametric}
