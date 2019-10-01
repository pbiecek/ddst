# Data Driven Smooth Tests with ddst package

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ddst)](http://cran.r-project.org/web/packages/ddst)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ddst?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/ddst)
[![Build Status](https://api.travis-ci.org/pbiecek/ddst.png)](https://travis-ci.org/pbiecek/ddst)

## Overview

`ddst` stands for `data driven smooth test`. The test characterizes data-dependent choice of the number of components in a smooth test statistic.

In this package you will find two groups of selected data driven smooth tests: goodness-of-fit tests and nonparametric tests for comparing distributions.
    
### Data Driven Smooth Tests for Selected Goodness-of-Fit Problems

These tests were inspired by the results from: [*Data driven smooth tests for composite hypotheses* by Inglot, Kallenberg, and Ledwina (1997)](https://projecteuclid.org/euclid.aos/1069362746) and [*Towards data driven selection of a penalty function for data driven Neyman tests* by Inglot and Ledwina (2006)](https://www.sciencedirect.com/science/article/pii/S0024379505005082).


* DDSTest for Uniformity - `ddst.uniform.test()`; see [*Towards data driven selection of a penalty function for data driven Neyman tests* by Inglot and Ledwina (2006)](https://www.sciencedirect.com/science/article/pii/S0024379505005082).
* DDSTest for Exponentiality - `ddst.exp.test()`; see [*Data driven smooth tests for composite hypotheses: Comparison of powers* by Kallenberg and Ledwina (1997)](https://www.tandfonline.com/doi/abs/10.1080/00949659708811850).
* DDSTest for Normality; bounded basis functions - `ddst.normbounded.test()`; see [*Data-driven tests for a location-scale family revisited* by Janic and Ledwina (2009)](https://link.springer.com/article/10.1080/15598608.2009.10411952).
* DDSTest for Normality; unbounded basis functions - `ddst.normunbounded.test()`; see [*Detection of non-Gaussianity* by Ledwina and Wylupek (2015)](https://www.tandfonline.com/doi/abs/10.1080/00949655.2014.983110?journalCode=gscs20).
* DDSTest for Extreme Value Distribution - `ddst.extr.test()`; see [*Data-driven tests for a location-scale family revisited* by Janic and Ledwina (2009)](https://link.springer.com/article/10.1080/15598608.2009.10411952).


### Nonparametric Data Driven Smooth Tests for Comparing Distributions

A starting point of the constructions were the papers: [*Data driven rank test for two-sample problem* by Janic-Wroblewska and Ledwina (2000)](https://www.jstor.org/stable/4616603?seq=1#page_scan_tab_contents) and [*Towards data driven selection of a penalty function for data driven Neyman tests* by Inglot and Ledwina (2006)](https://www.sciencedirect.com/science/article/pii/S0024379505005082).

* DDSTest for k-Sample Problem - `ddst.ksample.test()`; see [*Data-driven k-sample tests* by Wylupek (2010)](https://www.jstor.org/stable/40586684?seq=1).
* DDSTest for Change-Point Problem - `ddst.changepoint.test()`; see [*Data driven rank test for the change point problem* by Antoch, Huskova, Janic and Ledwina (2008)](https://link.springer.com/article/10.1007/s00184-007-0139-2).
* DDSTest for Stochastic Dominance in Two Samples - `ddst.forstochdom.test()`; see [*Nonparametric tests for stochastic ordering* by Ledwina and Wylupek (2012) ](https://link.springer.com/article/10.1007/s11749-011-0278-7).
* DDSTest Against Stochastic Dominance - `ddst.againststochdom.test()`; see [*Two-sample test against one-sided alternatives* by Ledwina and Wylupek (2012)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9469.2011.00787.x).
* DDSTest for Upward Trend - `ddst.upwardtrend.test()`; see [*Data-driven tests for trend* by Wylupek (2013)](https://www.tandfonline.com/doi/abs/10.1080/03610926.2012.697967).
* DDSTest for Umbrella Alternatives - `ddst.umbrella.test()`; see [*An automatic test for the umbrella alternatives* by Wylupek (2016)](https://onlinelibrary.wiley.com/doi/abs/10.1111/sjos.12231).

A more detailed overview is contained in [Data Driven Smooth Test - Introductory Material](http://www.biecek.pl/R/ddst/description.pdf). Full details on the above procedures can be found in the related papers.

## Installation

```{r}
# the easiest way to get ddst is to install it from CRAN:
install.packages("ddst")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("pbiecek/ddst")
```
