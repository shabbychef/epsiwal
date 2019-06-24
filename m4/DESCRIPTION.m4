dnl divert here just means the output from basedefs does not appear.
divert(-1)
include(basedefs.m4)
divert(0)dnl
Package: PKG_NAME()
Maintainer: Steven E. Pav <shabbychef@gmail.com>
Authors@R: c(person(c("Steven", "E."), "Pav", 
    role=c("aut","cre"),
    email="shabbychef@gmail.com",
    comment = c(ORCID = "0000-0002-4197-6195")))
Version: VERSION()
Date: DATE()
License: LGPL-3
Title: Exact Post Selection Inference with Applications to the Lasso
BugReports: https://github.com/shabbychef/PKG_NAME()/issues
Description: Implements the conditional estimation procedure of
  Lee, Sun, Sun and Taylor, <doi:10.1214/15-AOS1371>.
  This procedure allows hypothesis testing on the mean of
  a normal random vector subject to linear constraints.
Depends: 
    R (>= 3.0.2)
dnl Imports: 
dnl truncnorm
Suggests: 
dnl knitr,
    testthat
URL: https://github.com/shabbychef/PKG_NAME()
dnl VignetteBuilder: knitr
Collate:
m4_R_FILES()
dnl vim:ts=2:sw=2:tw=79:syn=m4:ft=m4:et
