% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ddst.evd.test.R
\name{ddst.evd.test}
\alias{ddst.evd.test}
\alias{ddst.extr.Nk}
\title{Data Driven Smooth Test for Extreme Value Distribution}
\usage{
ddst.evd.test(
  x,
  base = ddst.base.legendre,
  d.n = 10,
  c = 100,
  nr = 1e+05,
  compute.p = TRUE,
  alpha = 0.05,
  compute.cv = TRUE,
  ...
)
}
\arguments{
\item{x}{a (non-empty) numeric vector of data values}

\item{base}{a function which returns an orthonormal system, possible choice: \code{ddst.base.legendre} for the Legendre polynomials and \code{ddst.base.cos} for the cosine system}

\item{d.n}{an integer specifying the maximum dimension considered, only for advanced users}

\item{c}{a calibrating parameter in the penalty in the model selection rule}

\item{nr}{an integer specifying the number of runs for a p-value and a critical value computation if any}

\item{compute.p}{a logical value indicating whether to compute a p-value or not}

\item{alpha}{a significance level}

\item{compute.cv}{a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not}

\item{...}{further arguments}
}
\value{
An object of class \code{htest}
  \item{statistic }{the value of the test statistic.}
  \item{parameter }{the number of choosen coordinates (k).}
  \item{method }{a character string indicating the parameters of performed test. }
  \item{data.name }{a character string giving the name(s) of the data. }
  \item{p.value }{the p-value for the test, computed only if \code{compute.p=TRUE}.}
}
\description{
Performs data driven smooth test for composite hypothesis of extreme value distribution.
Null density is given by
\eqn{
  f(z;\gamma)=1/\gamma_2 \exp((z-\gamma_1)/\gamma_2- \exp((z-\gamma_1)/\gamma_2))}, \eqn{z \in R}.
}
\details{
We model alternatives similarly as in Kallenberg and Ledwina (1997) and Janic-Wroblewska (2004) using Legendre's polynomials or cosines.
  For more details see: \url{http://www.biecek.pl/R/ddst/description.pdf}.
}
\examples{
library(evd)
set.seed(7)

# for given vector of 19 numbers
z <- c(13.41, 6.04, 1.26, 3.67, -4.54, 2.92, 0.44, 12.93, 6.77, 10.09,
      4.10, 4.04, -1.97, 2.17, -5.38, -7.30, 4.75, 5.63, 8.84)
\dontrun{
t <- ddst.evd.test(z, compute.p = TRUE, d.n = 10)
t
plot(t)

# H0 is true
x <- -qgumbel(runif(100),-1,1)
t <- ddst.evd.test (x, compute.p = TRUE, d.n = 10)
t
plot(t)

# H0 is false
x <- rexp(80,4)
t <- ddst.evd.test (x, compute.p = TRUE, d.n = 10)
t
plot(t)
}
}
\references{
Hosking, J.R.M., Wallis, J.R., Wood, E.F. (1985). Estimation of the generalized extreme-value distribution by the method of probability-weighted moments. \eqn{ Technometrics} 27, 251--261.

Janic-Wroblewska, A. (2004). Data-driven smooth test for extreme  value distribution. \eqn{ Statistics} 38, 413--426.

Janic, A. and Ledwina, T. (2008). Data-driven tests for a location-scale family revisited. \eqn{ J. Statist. Theory. Pract. Special issue on Modern Goodness of Fit Methods. accepted.}.

Kallenberg, W.C.M., Ledwina, T. (1997). Data driven smooth tests for composite hypotheses: Comparison of powers. \eqn{ J. Statist. Comput. Simul.} \bold{ 59}, 101--121.
}
\keyword{htest}
