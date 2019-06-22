

# epsiwal

[![Build Status](https://travis-ci.org/shabbychef/epsiwal.png)](https://travis-ci.org/shabbychef/epsiwal)
[![codecov.io](http://codecov.io/github/shabbychef/epsiwal/coverage.svg?branch=master)](http://codecov.io/github/shabbychef/epsiwal?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/epsiwal)](https://cran.r-project.org/package=epsiwal)
[![Downloads](http://cranlogs.r-pkg.org/badges/epsiwal?color=green)](http://www.r-pkg.org/pkg/epsiwal)
[![Total](http://cranlogs.r-pkg.org/badges/grand-total/epsiwal?color=green)](http://www.r-pkg.org/pkg/epsiwal)

Implements the conditional inference on normal variates described in 
Lee, Sun, Sun and Taylor, "Exact Post Selection Inference, with Application to the Lasso."


-- Steven E. Pav, shabbychef@gmail.com

## Installation

This package may be installed from CRAN; the latest version may be
found on [github](https://www.github.com/shabbychef/epsiwal "epsiwal")
via devtools, or installed via [drat](https://github.com/eddelbuettel/drat "drat"):


```r
# CRAN
install.packages(c("epsiwal"))
# devtools
if (require(devtools)) {
    # latest greatest
    install_github("shabbychef/epsiwal")
}
# via drat:
if (require(drat)) {
    drat:::add("shabbychef")
    install.packages("epsiwal")
}
```

# Basic Usage


```r
requre(epsiwal)
```

