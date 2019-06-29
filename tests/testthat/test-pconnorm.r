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
# Created: 2019.06.21
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) { set.seed(as.integer(charToRaw(str))) }
#UNFOLD

context("pconnorm")
test_that("runs without error",{#FOLDUP
	set.char.seed("d106bf01-7ad5-49b2-8715-bb8148a0e294")

	set.seed(1234)
	n <- 10
	y <- rnorm(n)
	A <- matrix(rnorm(n*(n-3)),ncol=n)
	b <- A%*%y + runif(nrow(A))
	Sigma <- diag(runif(n))
	mu <- rnorm(n)
	eta <- rnorm(n)

	# can we use it?
	expect_error(pval <- pconnorm(y=y,A=A,b=b,eta=eta,mu=mu,Sigma=Sigma),NA)

	# lower.tail
	expect_error(mpval <- pconnorm(y=y,A=A,b=b,eta=eta,mu=mu,Sigma=Sigma,lower.tail=FALSE),NA)
	expect_equal(mpval,1-pval,tolerance=1e-5)

	# log p
	expect_error(lval <- pconnorm(y=y,A=A,b=b,eta=eta,mu=mu,Sigma=Sigma,log.p=TRUE),NA)
	expect_equal(lval,log(pval),tolerance=1e-5)

	# premultiply
	expect_error(pval1 <- pconnorm(y=y,A=A,b=b,eta=eta,mu=mu,Sigma=Sigma),NA)
	expect_error(pval2 <- pconnorm(y=y,A=A,b=b,eta=eta,mu=mu,Sigma_eta=Sigma %*% eta),NA)
	expect_error(pval3 <- pconnorm(y=y,A=A,b=b,eta=eta,Sigma=Sigma,eta_mu=as.numeric(t(eta)%*%mu)),NA)
	expect_error(pval4 <- pconnorm(y=y,A=A,b=b,eta=eta,Sigma_eta=Sigma %*% eta,eta_mu=as.numeric(t(eta)%*%mu)),NA)
	expect_equal(pval1,pval2,tolerance=1e-5)
	expect_equal(pval1,pval3,tolerance=1e-5)
	expect_equal(pval1,pval4,tolerance=1e-5)

	# pval is decreasing in eta mu 
	aleta <- y
	expect_error(pval_lo <- pconnorm(y=y,A=A,b=b,eta=aleta,Sigma_eta=Sigma%*%aleta,eta_mu=as.numeric(t(aleta)%*%mu)),NA)
	expect_error(pval_hi <- pconnorm(y=y,A=A,b=b,eta=aleta,Sigma_eta=Sigma%*%aleta,eta_mu=as.numeric(t(aleta)%*%mu)-1),NA)
	expect_lt(pval_lo,pval_hi)
})#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
