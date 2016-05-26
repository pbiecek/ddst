\name{ddst.extr.test}
\alias{ddst.extr.test}
\alias{ddst.extr.Nk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Data Driven Smooth Test for Extreme Value Distribution }
\description{
Performs data driven smooth test for composite hypothesis of extreme value distribution.
}
\usage{
ddst.extr.test(x, base = ddst.base.legendre, c = 100, B = 1000, compute.p = F, 
    Dmax = 5, ...)
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
Null density is given by 
\emph{
$f(z;gamma)=1/gamma_2 exp((z-gamma_1)/gamma_2- exp((z-gamma_1)/gamma_2))$}, z in R.

We model alternatives similarly as in Kallenberg and Ledwina (1997) and Janic-Wroblewska (2004) using Legendre's polynomials or cosines. The parameter 
\emph{$gamma=(gamma_1,gamma_2)$} is estimated by \emph{$tilde gamma=(tilde gamma_1,tilde gamma_2)$}, where \emph{$tilde gamma_1=-1/n sum_{i=1}^n Z_i + varepsilon  G$}, where \emph{$varepsilon approx 0.577216 $} is the Euler constant and \emph{$ G = tilde gamma_2 = [n(n-1) ln2]^{-1}sum_{1<= j < i <= n}(Z_{n:i}^o - Z_{n:j}^o) $} while \emph{$Z_{n:1}^o <= ... <= Z_{n:n}^o$}
are ordered variables \emph{$-Z_1,...,-Z_n$}, cf Hosking et al. (1985). 
The above yields auxiliary test statistic \emph{$W_k^*(tilde gamma)$} described in details in Janic and Ledwina (2008), in case when Legendre's basis is applied. 

The related matrix \emph{$[I^*(tilde gamma)]^{-1}$} does not  depend on \emph{$tilde gamma$} and is calculated for succeding dimensions \emph{k} using some recurrent relations for Legendre's polynomials and numerical methods for cosine functions. In the implementation the default value of \emph{c} in \emph{$T^*$} was fixed to be 100. Hence, \emph{$T^*$} is Schwarz-type model selection rule. The resulting data driven test statistic for extreme value distribution is \emph{$W_{T^*}=W_{T^*}(tilde gamma)$}.

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
Hosking, J.R.M., Wallis, J.R., Wood, E.F. (1985). Estimation of the generalized extreme-value distribution by the method of probability-weighted moments. \emph{ Technometrics} 27, 251--261.

Janic-Wroblewska, A. (2004). Data-driven smooth test for extreme  value distribution. \emph{ Statistics} 38, 413--426.

Janic, A. and Ledwina, T. (2008). Data-driven tests for a location-scale family revisited. \emph{ J. Statist. Theory. Pract. Special issue on Modern Goodness of Fit Methods. accepted.}.

Kallenberg, W.C.M., Ledwina, T. (1997). Data driven smooth tests for composite hypotheses: Comparison of powers. \emph{ J. Statist. Comput. Simul.} \bold{ 59}, 101--121. 
}
\author{ Przemyslaw Biecek and Teresa Ledwina }

\examples{
library(evd)

# for given vector of 19 numbers
z = c(13.41, 6.04, 1.26, 3.67, -4.54, 2.92, 0.44, 12.93, 6.77, 10.09, 
   4.10, 4.04, -1.97, 2.17, -5.38, -7.30, 4.75, 5.63, 8.84)
ddst.extr.test(z, compute.p=TRUE)

# H0 is true
x = -qgumbel(runif(100),-1,1)
ddst.extr.test (x, compute.p = TRUE)

# H0 is false
x = rexp(80,4)
ddst.extr.test (x, compute.p = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{nonparametric}
