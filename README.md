Data Driven Smooth Tests with ddst package
==========================================

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ddst)](http://cran.r-project.org/web/packages/ddst)
[![Downloads](http://cranlogs.r-pkg.org/badges/ddst)](http://cran.rstudio.com/package=ddst)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ddst?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/ddst)

`ddst` stands for Data Driven Smooth Tests, these are tests in which alternative hypothesis is dynamically selected based on data.

In this package you will find two groups of smooth of test: goodness-of-fit tests and
    nonparametric tests for comparing distributions.
    
Among goodness-of-fit tests there are tests for

* `ddst.exp.test` DDSTest for Exponentiality, find more: http://pbiecek.github.io/ddst/reference/ddst.exp.test.html, 
* `ddst.extr.test` DDSTest for  Extreme Value Distribution, find more: http://pbiecek.github.io/ddst/reference/ddst.extr.test.html, 
* `ddst.norm.test` DDSTest for Normality, find more: http://pbiecek.github.io/ddst/reference/ddst.norm.test.html, 
* `ddst.uniform.test` DDSTest for Uniformity, find more: http://pbiecek.github.io/ddst/reference/ddst.uniform.test.html, 

Among nonparametric tests there are tests for 

* `ddst.ksample.test` DDSTest for k-Sample Test, find more: http://pbiecek.github.io/ddst/reference/ddst.ksample.test.html, 
* `ddst.umbrella.test` DDSTest for k-Sample Umbrella Alternatives, find more: http://pbiecek.github.io/ddst/reference/ddst.umbrella.test.html, 
* `ddst.twosample.test` DDSTest for Two-Sample Test Against One-Sided Alternatives, find more: http://pbiecek.github.io/ddst/reference/ddst.twosample.test.html, 
* `ddst.changepoint.test` DDSTest for Change-point problems,  
* `ddst.dominance.test` DDSTest for Stochastic dominance 

Find a detailed description for methofdology behond smooth tests at http://www.biecek.pl/R/ddst/description.pdf (TODO: update).

# How to install the ddst package

Install latest stable version of the `ddst` package from CRAN with

```{r}
install.packages("ddst")
```

or install development version version from GitHub 

```{r}
devtools::install_github("pbiecek/ddst")
```

