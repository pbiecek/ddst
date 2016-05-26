\name{ddst.exp.test}
\alias{ddst.exp.test}
\alias{ddst.exp.Nk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Data Driven Smooth Test for Exponentiality}
\description{
Performs data driven smooth test for composite hypothesis of exponentiality.
}
\usage{
ddst.exp.test(x, base = ddst.base.legendre, c = 100, B = 1000, compute.p = F, 
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

Null density is given by \emph{$f(z;gamma) = exp(-z/gamma)$}  for z >= 0 and 0 otherwise.  

Modelling alternatives similarly as in Kallenberg and Ledwina (1997 a,b), e.g., and estimating \emph{$gamma$} by \emph{$tilde gamma= 1/n sum_{i=1}^n Z_i$} yields the efficient score 
vector \emph{$l^*(Z_i;tilde gamma)=(phi_1(F(Z_i;tilde gamma)),...,phi_k(F(Z_i;tilde gamma)))$}, where \emph{$phi_j$}'s are	\emph{j}th degree orthonormal Legendre polynomials on [0,1] or cosine functions \emph{$sqrt(2) cos(pi j x), j>=1$,} while \emph{$F(z;gamma)$} is the distribution function pertaining to \emph{$f(z;gamma)$}. 

The matrix \emph{$[I^*(tilde gamma)]^{-1}$} does not  depend on \emph{$tilde gamma$} and is calculated for succeding dimensions \emph{k} using some recurrent relations for Legendre's polynomials and computed in a numerical way in case of cosine basis. In the implementation the default value of \emph{c} in \emph{$T^*$} is set to be 100. 

Therefore, \emph{$T^*$} practically coincides with S1 considered in Kallenberg and Ledwina (1997 a).

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
Kallenberg, W.C.M., Ledwina, T. (1997 a). Data driven smooth tests for
composite hypotheses: Comparison of powers. \emph{ J. Statist. Comput. Simul.}
\bold{ 59}, 101--121.

Kallenberg, W.C.M.,  Ledwina, T. (1997 b). Data driven smooth tests when
the hypothesis is composite. \emph{ J. Amer. Statist. Assoc.} \bold{ 92},
1094--1104.

}
\author{ Przemyslaw Biecek and Teresa Ledwina }
\examples{

# H0 is true
z = rexp(80,4)
ddst.exp.test (z, compute.p = TRUE)

# H0 is false
z = rchisq(80,4)
(t = ddst.exp.test (z, compute.p = TRUE))
t$p.value

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{nonparametric}
