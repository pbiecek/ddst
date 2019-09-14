# Data Driven Smooth Tests with ddst package

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ddst)](http://cran.r-project.org/web/packages/ddst)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ddst?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/ddst)
[![Build Status](https://api.travis-ci.org/pbiecek/ddst.png)](https://travis-ci.org/pbiecek/ddst)

## Overview

`ddst` stands for Data Driven Smooth Tests. These tests characterize data-dependent choice of the number of component in test statistics.

In this package you will find two groups of data driven smooth tests: goodness-of-fit tests and nonparametric tests for comparing distributions.
    
### Data Driven Smooth Tests for Goodness-of-fit Problems

These tests are based mostly on [Data-driven version of Neyman's smooth test of fit. Ledwina (1994)](https://www.jstor.org/stable/2290926?seq=1).

* DDSTest for Uniformity - `ddst.uniform.test()`.
* DDSTest for Exponentiality - `ddst.exp.test()`.
* DDSTest for Extreme Value Distribution - `ddst.extr.test()`.
* DDSTest for Normality for Bounded basis functions - `ddst.norm.test()`.
* DDSTest for Normality, Unbounded basis functions - `ddst.normunbounded.test()`. This test is described in details in [Detection of non-Gaussianity. Ledwina and Wylupek (2015)](https://www.tandfonline.com/doi/abs/10.1080/00949655.2014.983110?journalCode=gscs20).

### Nonparametric Data Driven Smooth Tests for Comparing Distributions

* DDSTest for k-sample problem - `ddst.ksample.test()`. This test is described in details in [Data-driven k-sample tests. Wylupek (2010)](https://www.jstor.org/stable/40586684?seq=1).
* DDSTest for two-sample problem - `ddst.twosample.test()`. This test is a specific version of the k-sample test.
* DDSTest for nonparametric change point problem - `ddst.changepoint.test()`. This test is described in details in [Data driven rank test for the change point problem. Antoch, Huskova, Janic and Ledwina (2008)](https://link.springer.com/article/10.1007/s00184-007-0139-2).
* DDSTest for detection of stochastic dominance of the first order - `ddst.stochasticorder.test()`. This test is described in details in [Nonparametric tests for stochastic ordering. Ledwina and Wylupek (2012) ](https://link.springer.com/article/10.1007/s11749-011-0278-7).
* DDSTest for stochastic dominance - `ddst.stochasticdominance.test()`. This test is described in details in [Two-sample test against one-sided alternatives. Ledwina and Wylupek (2012)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9469.2011.00787.x).
* DDSTest for the two- and k-sample problems with stochastically ordered alternatives - `ddst.upward.test()`. This test is described in details in [Data-driven tests for trend. Wylupek (2013)](https://www.tandfonline.com/doi/abs/10.1080/03610926.2012.697967).
* DDSTest for k-sample umbrella alternatives - `ddst.umbrella.test()`. This test is described in details in [An automatic test for the umbrella alternatives. Wylupek (2016)](https://onlinelibrary.wiley.com/doi/abs/10.1111/sjos.12231).

## Installation

```{r}
# the easiest way to get ddst is to install it from CRAN:
install.packages("ddst")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("pbiecek/ddst")
```
