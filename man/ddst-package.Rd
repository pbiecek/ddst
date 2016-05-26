\name{ddst-package}
\alias{ddst-package}
\alias{ddst.base.cos}
\alias{ddst.base.legendre}
\alias{ddst.IIC}
\alias{ddst.phi}
\alias{ddst}
\alias{Nmax}
\docType{package}
\title{
Data Driven Smooth Tests
}
\description{
Set of Data Driven Smooth Tests for Goodness of Fit
}
\details{
\tabular{ll}{
Package: \tab ddst\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2008-07-01\cr
License: \tab GPL-2\cr
}
}
\section{General Description}{
Smooth test was introduced by Neyman (1937) to verify simple null hypothesis asserting that observations obey completely known continuous distribution function \emph{F}. Smooth test statistic (with \emph{k} components) can be interpreted as score statistic in an appropriate class of auxiliary models indexed by a vector of parameters \emph{$theta in R^k, k >= 1$.} 

Pertaining auxilary null hypothesis asserts \emph{$theta=theta_0=0$}. Therefore, in this case, the smooth test statistic based on \emph{n} i.i.d. observations \emph{$Z_1,...,Z_n$} has the form
\emph{$W_k=[1/sqrt(n) sum_{i=1}^n l(Z_i)]I^{-1}[1/sqrt(n) sum_{i=1}^n l(Z_i)]'$},

where \emph{$l(Z_i)$}, i=1,...,n, is \emph{k}-dimensional (row) score vector, the symbol \emph{'} denotes transposition while \emph{$I=Cov_{theta_0}[l(Z_1)]'[l(Z_1)]$}. Following Neyman's idea of modelling underlying distributions one gets \emph{$l(Z_i)=(phi_1(F(Z_i)),...,phi_k(F(Z_i)))$} and \emph{I} being the identity matrix, where \emph{$phi_j$}'s, j >= 1,  are zero mean orthonormal functions on [0,1], while \emph{F} is the completely specified null distribution function.

In case of composite null hypothesis there is also unspecified vector of nuisance parameters \emph{$gamma$} defining the distribution of observations. Smooth statistic (with \emph{k} components) in such applications is understood as efficient score statistic for some class of models indexed by an auxiliary parmeter \emph{$theta in R^k$}, k >= 1. Pertaining efficient score vector \emph{$l^*(Z_i;gamma)$} is defined as the residual from projection the score vector for \emph{$theta$} onto the space spanned by score vector for \emph{$gamma$}. As such, smooth test is alternative name for \emph{$C(alpha)$} Neyman's test. See Neyman (1959), Buhler and Puri (1966) as well as Javitz (1975) for details. Hence, smooth test, based on \emph{n} i.i.d. variables \emph{$Z_1,...,Z_n$} rejects hypothesis \emph{$theta=theta_0=0$} for large values of 

\emph{$W_k^{*}(tilde gamma)=[1/sqrt(n) sum_{i=1}^n l^*(Z_i;tilde gamma)][I^*(tilde gamma)]^{-1}[1/sqrt(n) sum_{i=1}^n l^*(Z_i;tilde gamma)]'$},
where \emph{$tilde gamma$} is an appropriate estimator of \emph{$gamma$} while \emph{$I^*(gamma)=Cov_{theta_0}[l^*(Z_1;gamma)]'[l^*(Z_1;gamma)]$}. More details can be found in Janic and Ledwina (2008), Kallenberg and Ledwina (1997 a,b) as well as Inglot and Ledwina (2006 a,b). 

Auxiliary models, mentioned above, aim to mimic the unknown underlying model for the data at hand. To choose the dimension \emph{k} of the auxilary model we apply some model selection criteria. Among several solutions already considered, we decided to implement two following ones, pertaining to the two above described problems and resulting \emph{$W_k$} and \emph{$W_k^*(tilde gamma)$}. The selection rules in the two cases are briefly denoted by \emph{T} and \emph{$T^*$}, respectively, and given by

\emph{$T = min{1 <= k <= d: W_k-pi(k,n,c) >= W_j-pi(j,n,c), j=1,...,d}$}

and

\emph{
$T^* = min{1 <= k <= d: W_k^*(tilde gamma)-pi^*(k,n,c) >= W_j^*(tilde gamma)-pi^*(j,n,c), j=1,...,d}$}.

Both criteria are based on approximations of penalized loglikelihoods, where loglikelihoods are replaced by \emph{$W_k$} and \emph{$W_k^*(tilde gamma)$}, respectively.  The penalties for the dimension \emph{j} in case of simple and composite null hypothesis are defined as follows

\emph{$pi(j,n,c)={jlog n,  if  max{1 <= k <= d}|Y_k| <= sqrt(c log(n)), 2j,  if max{1 <= k <= d}|Y_k|>sqrt(c log(n)). }$}

and

\emph{
$pi^*(j,n,c)={jlog n,  if max{1 <= k <= d}|Y_k^*| <= sqrt(c log(n)),2j  if max(1 <= k <= d)|Y_k^*| > sqrt(c log(n))}$}.

respectively, where \emph{c} is some calibrating constant, \emph{d} is maximal dimension taken into account, 

\emph{$(Y_1,...,Y_k)=[1/sqrt(n) sum_{i=1}^n l(Z_i)]I^{-1/2}$} 

while  

\emph{$(Y_1^*,...,Y_k^*)=[1/sqrt(n) sum_{i=1}^n l^*(Z_i; tilde gamma)][I^*(tilde gamma)]^{-1/2}$}.

In consequence, data driven smooth tests for the simple and composite null hypothesis reject for large values of \emph{$W_T$}
and \emph{$W_{T^*} = W_{T^*}(tilde gamma)$}, respectively. For details see Inglot and Ledwina (2006 a,b,c). 

The choice of \emph{c} in \emph{T} and \emph{$T^*$} is decisive to finite sample behaviour of the selection rules and pertaining statistics \emph{$W_T$} and \emph{$W_{T^*}(tilde gamma)$}. In particular, under large \emph{c}'s the rules behave similarly as Schwarz's (1978) BIC while for \emph{c=0} they mimic Akaike's (1973) AIC. For moderate sample sizes, values \emph{c in (2,2.5)} guarantee, under `smooth' departures, only slightly smaller power as in case BIC were used and simultaneously give much higher power than BIC under multimodal alternatives. In genral, large \emph{c's} are recommended if changes in location, scale, skewness and kurtosis are in principle aimed to be detected. For evidence and discussion see Inglot and Ledwina (2006 c). 

It \emph{c>0} then the limiting null distribution of \emph{$W_T$} and \emph{$W_{T^*}(tilde gamma)$} is central chi-squared with one degree of freedom. In our implementation, for given \emph{n}, both critical values and \emph{p}-values are computed by MC method.
 
Empirical distributions of \emph{T} and \emph{$T^*$} as well as \emph{$W_T$} and \emph{$W_{T^*}(tilde gamma)$} are not essentially influenced by the choice of reasonably large \emph{d}'s, provided that sample size is at least moderate.

For more details see: \url{http://www.biecek.pl/R/ddst/description.pdf}.

}
\author{
Przemyslaw Biecek and Teresa Ledwina

Maintainer: You should complain to Przemyslaw Biecek  <przemyslaw.biecek@gmail.com>
}
\references{
Akaike, H. (1973). Information theory and the maximum likelihood principle. In: \emph{ 2nd International Symposium on Information Theory}, (eds. B. N. Petrov and F. Csaki), 267-281. Akademiai Kiado, Budapest.

Buhler, W.J., Puri, P.S. (1966). On optimal asymptotic tests of composite hypotheses with several constraints. \emph{ Z. Wahrsch. verw. Geb.} \bold{ 5}, 71--88.

Inglot, T., Ledwina, T. (2006 a). Data-driven score tests for homoscedastic linear regression model: asymptotic results. \emph{ Probab. Math. Statist.} \bold{ 26}, 41--61.

Inglot, T., Ledwina, T. (2006 b). Data-driven score tests for homoscedastic linear regression model: the construction and simulations. In \emph{ Prague Stochastics 2006. Proceedings}, (eds. M. Huskova, M. Janzura), 124--137. Matfyzpress, Prague.

Inglot, T., Ledwina, T. (2006 c). Towards data driven selection of a penalty function for data driven Neyman tests. \emph{ Linear Algebra and its Appl.} \bold{ 417}, 579--590. 

Javitz, H.S. (1975). Generalized smooth tests of goodness of fit, independence and equality of distributions. Ph.D. thesis at University of California, Berkeley.

Janic, A. and Ledwina, T. (2008). Data-driven tests for a location-scale family revisited. \emph{ J. Statist. Theory. Pract. Special issue on Modern Goodness of Fit Methods. accepted.}.

Kallenberg, W.C.M., Ledwina, T. (1997 a). Data driven smooth tests for composite hypotheses: Comparison of powers. \emph{ J. Statist. Comput. Simul.} \bold{ 59}, 101--121.

Kallenberg, W.C.M.,  Ledwina, T. (1997 b). Data driven smooth tests when the hypothesis is composite. \emph{ J. Amer. Statist. Assoc.} \bold{ 92}, 1094--1104.

Neyman, J. (1937). `Smooth test' for goodness of fit. \emph{ Skand. Aktuarietidskr.} \bold{ 20}, 149-199.

Neyman, J. (1959). Optimal asymptotic tests of composite statistical hypotheses. In \emph{ Probability and Statistics}, (ed. U. Grenander), Harald Cramer Volume, 212--234. Wiley, New York.
}
\keyword{ package }
\keyword{htest}
\keyword{nonparametric}
\examples{

# Data Driven Smooth Test for Uniformity
#
# H0 is true
z = runif(80)
ddst.uniform.test(z, compute.p=TRUE)

# H0 is false
z = rbeta(80,4,2)
(t = ddst.uniform.test(z, compute.p=TRUE))
t$p.value

# Data Driven Smooth Test for Normality
#
# H0 is true
z = rnorm(80)
ddst.norm.test(z, compute.p=TRUE)

# H0 is false
z = rexp(80,4)
ddst.norm.test(z, B=5000, compute.p=TRUE)

# Data Driven Smooth Test for Extreme Value Distribution
#
# H0 is true
#library(evd)
#z = -qgumbel(runif(100),-1,1)
#ddst.extr.test (z, compute.p = TRUE)

# H0 is false
z = rexp(80,4)
ddst.extr.test (z, compute.p = TRUE)

# Data Driven Smooth Test for Exponentiality
#
# H0 is true
z = rexp(80,4)
ddst.exp.test (z, compute.p = TRUE)

# H0 is false
z = rchisq(80,4)
ddst.exp.test (z, compute.p = TRUE)

}
