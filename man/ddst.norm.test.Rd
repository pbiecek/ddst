\name{ddst.norm.test}
\alias{ddst.norm.test}
\alias{tabNorm}
\alias{ddst.norm.Nk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Data Driven Smooth Test for Normality }
\description{
Performs data driven smooth test for composite hypothesis of normality.
}
\usage{
ddst.norm.test(x, base = ddst.base.legendre, c = 100, B = 1000, compute.p = F, 
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
$f(z;gamma)=1/(sqrt(2 pi)gamma_2) exp(-(z-gamma_1)^2/(2 gamma_2^2))$} for z in R.

We model alternatives similarly as in Kallenberg and Ledwina (1997 a,b) using Legendre's polynomials or cosine basis. The parameter 
\emph{$gamma=(gamma_1,gamma_2)$} is estimated by \emph{$tilde gamma=(tilde gamma_1,tilde gamma_2)$}, where \emph{$tilde gamma_1=1/n sum_{i=1}^n Z_i$} and 
\emph{$tilde gamma_2 = 1/(n-1) sum_{i=1}^{n-1}(Z_{n:i+1}-Z_{n:i})(H_{i+1}-H_i)$},
while \emph{$Z_{n:1}<= ... <= Z_{n:n}$} are ordered values of \emph{$Z_1, ..., Z_n$} and \emph{$H_i= phi^{-1}((i-3/8)(n+1/4))$}, cf. Chen and Shapiro (1995). 

The above yields auxiliary test statistic \emph{$W_k^*(tilde gamma)$} described in details in Janic and Ledwina (2008), in case when Legendre's basis is applied. 
The pertaining matrix \emph{$[I^*(tilde gamma)]^{-1}$} does not  depend on \emph{$tilde gamma$} and is calculated for succeding dimensions \emph{k} using some recurrent relations for Legendre's polynomials and is computed in a numerical way in case of cosine basis. In the implementation of \emph{$T^*$} the default value of \emph{c} is set  to be 100. Therefore, in practice, \emph{$T^*$} is Schwarz-type criterion. See Inglot and Ledwina (2006) as well as Janic and Ledwina (2008) for comments. The resulting data driven test statistic for normality is \emph{$W_{T^*}=W_{T^*}(tilde gamma)$}.

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
Chen, L., Shapiro, S.S. (1995). An alternative test for normality based on normalized spacings. \emph{ J. Statist. Comput. Simulation} 53, 269--288.

Inglot, T., Ledwina, T. (2006). Towards data driven selection of a penalty function for data driven Neyman tests. \emph{ Linear Algebra and its Appl.} \bold{ 417}, 579--590. 

Janic, A. and Ledwina, T. (2008). Data-driven tests for a location-scale family revisited. \emph{ J. Statist. Theory. Pract. Special issue on Modern Goodness of Fit Methods. accepted.}.

Kallenberg, W.C.M., Ledwina, T. (1997 a). Data driven smooth tests for composite hypotheses: Comparison of powers. \emph{ J. Statist. Comput. Simul.} \bold{ 59}, 101--121.

Kallenberg, W.C.M.,  Ledwina, T. (1997 b). Data driven smooth tests when the hypothesis is composite. \emph{ J. Amer. Statist. Assoc.} \bold{ 92}, 1094--1104.

 }
\author{ Przemyslaw Biecek and Teresa Ledwina }
\examples{

# for given vector of 19 numbers
z = c(13.41, 6.04, 1.26, 3.67, -4.54, 2.92, 0.44, 12.93, 6.77, 10.09, 
   4.10, 4.04, -1.97, 2.17, -5.38, -7.30, 4.75, 5.63, 8.84)
ddst.norm.test(z, compute.p=TRUE)

# H0 is true
z = rnorm(80)
ddst.norm.test(z, compute.p=TRUE)

# H0 is false
z = rexp(80,4)
ddst.norm.test(z, B=5000, compute.p=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{nonparametric}
