# /usr/bin/r
#
# Copyright 2019-2019 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav 
#
# This file is part of epsiwal.
#
# epsiwal is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# epsiwal is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with epsiwal.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2019.06.23
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) { set.seed(as.integer(charToRaw(str))) }
#UNFOLD

context("ci_pconnorm")
test_that("runs without error",{#FOLDUP
	set.char.seed("fc742d8c-3454-440a-afd6-0bfd9d27e4c5")

	set.seed(1234)
	n <- 10
	y <- rnorm(n)
	A <- matrix(rnorm(n*(n-3)),ncol=n)
	b <- A%*%y + runif(nrow(A))
	Sigma <- diag(runif(n))
	mu <- rnorm(n)
	eta <- rnorm(n)


	pvals <- seq(0,1,by=0.1)
	expect_error(civals <- ci_connorm(y=y,A=A,b=b,eta=eta,Sigma=Sigma,p=pvals),NA)
	expect_equal(length(pvals),length(civals))
	expect_true(all(diff(civals) < 0))

	# integration tests:
	expect_error(pval <- pconnorm(y=y,A=A,b=b,eta=eta,mu=mu,Sigma=Sigma),NA)
	expect_error(cival <- ci_connorm(y=y,A=A,b=b,eta=eta,Sigma=Sigma,p=pval),NA)
	expect_equal(cival,as.numeric(t(eta) %*% mu),tolerance=1e-4)
})#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
