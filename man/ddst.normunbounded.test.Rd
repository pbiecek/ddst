% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ddst.normunbounded.test.R
\name{ddst.normubounded.test}
\alias{ddst.normubounded.test}
\alias{ddst.normunbounded.bias}
\title{Data Driven Smooth Test for Normality; Unbounded Basis Functions}
\usage{
ddst.normubounded.test(
  x,
  d.n = 20,
  e.0,
  v.0,
  r.alpha,
  s.n.alpha,
  alpha = 0.05,
  nr = 10000,
  compute.cv = TRUE
)

ddst.normunbounded.bias(n = 100, d.n = 20, nr = 10000)
}
\arguments{
\item{x}{a (non-empty) numeric vector of data values}

\item{d.n}{an integer specifying the maximum dimension considered, only for advanced users}

\item{e.0}{a (non-empty) numeric vector being the Monte Carlo estimate of the mean of the vector (C2; ... ;Cd.n) calculated using the function \code{ddst.normunbounded.bias()$e.0}}

\item{v.0}{a (non-empty) numeric vector being the Monte Carlo estimate of the variance of the vector (C2; ... ;Cd.n) calculated using the function \code{ddst.normunbounded.bias()$e.0}}

\item{r.alpha}{a critical value of the alpha level R.n test}

\item{s.n.alpha}{a penalty in the auxiliary model selection rule}

\item{alpha}{a significance level}

\item{nr}{an integer specifying the number of runs for a critical value computation}

\item{compute.cv}{a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not}

\item{n}{sample size}
}
\description{
Performs data driven smooth test for composite hypothesis of normality.
Null density is given by
\eqn{
  f(z;\gamma)=1/(\sqrt{2 \pi}\gamma_2) \exp(-(z-\gamma_1)^2/(2 \gamma_2^2))} for \eqn{z \in R}.
}
\examples{
set.seed(7)
# H0 is true
z <- rnorm(100)
# let's look on first 20 coordinates
d.n <- 20

\dontrun{
# calculate finite sample corrections
# see 6.2. Composite null hypothesis H in the appendix materials
e.v <- ddst.normunbounded.bias(n = length(z))
e.v

# simulated 1-alpha qunatiles, s(n, alpha) and  s.o(n, alpha)
# see Table 1 in the JSCS article
s <- 4.4
r.alpha <- 2.708

t <- ddst.normubounded.test(z, d.n, e.v$e.0, e.v$v.0, r.alpha, s)
t
plot(t)

# for Tephra data
z <- c(-1.748789, -1.75753, -1.740102, -1.740102, -1.731467, -1.765523,
       -1.761521, -1.72522, -1.80371, -1.745624, -1.872957, -1.729121,
       -1.81529, -1.888637, -1.887761, -1.881645, -1.91518, -1.849769,
       -1.755141, -1.665687, -1.764721, -1.736171, -1.736956, -1.737742,
       -1.687537, -1.804534, -1.790593, -1.808661, -1.784081, -1.729903,
       -1.711263, -1.748789, -1.772755, -1.72756, -1.71358, -1.821116,
       -1.839588, -1.839588, -1.830321, -1.807835, -1.747206, -1.788147,
       -1.759923, -1.786519, -1.726779, -1.738528, -1.754345, -1.781646,
       -1.641949, -1.755936, -1.775175, -1.736956, -1.705103, -1.743255,
       -1.82613, -1.826967, -1.780025, -1.684504, -1.751168)

# calculate finite sample corrections
e.v <- ddst.normunbounded.bias(n = length(z))
e.v

s <- 3.3
so <- 3.6
r.alpha <- 2.142

t <- ddst.normubounded.test(z, d.n, e.v$e.0, e.v$v.0, r.alpha, s)
t
plot(t)
}
}
\references{
Ledwina, T., Wyłupek, G. (2015) Detection of non-Gaussianity by Ledwina and Wyłupek \emph{Journal of Statistical Computation and Simulation} 17, 3480-3497.
}
\keyword{htest}
