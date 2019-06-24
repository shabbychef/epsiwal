

# epsiwal

[![Build Status](https://travis-ci.org/shabbychef/epsiwal.png)](https://travis-ci.org/shabbychef/epsiwal)
[![codecov.io](http://codecov.io/github/shabbychef/epsiwal/coverage.svg?branch=master)](http://codecov.io/github/shabbychef/epsiwal?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/epsiwal)](https://cran.r-project.org/package=epsiwal)
[![Downloads](http://cranlogs.r-pkg.org/badges/epsiwal?color=green)](http://www.r-pkg.org/pkg/epsiwal)
[![Total](http://cranlogs.r-pkg.org/badges/grand-total/epsiwal?color=green)](http://www.r-pkg.org/pkg/epsiwal)

Implements conditional inference on normal variates as described in 
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
library(epsiwal)
```

## See also

	* The original paper, by Lee, J. D., Sun, D. L., Sun, Y. and Taylor, J. E. 
		[Exact post-selection inference, with application to the Lasso](https://arxiv.org/abs/1311.6238).
	* The [`PSAT` package](https://github.com/ammeir2/PSAT), which supports similar
		procedures, but is not yet on CRAN.
	* The [`SelectiveInference` package](https://cran.r-project.org/web/packages/selectiveInference/index.html),
		which implements similar inferential procedures under quadratic
		constraints, as detailed in Tibshirani, R. J., Taylor, J., Lockhart, R. and Tibshirani, R. 
		[Exact Post-Selection Inference for Sequential Regression Procedures](https://arxiv.org/abs/1401.3889).

