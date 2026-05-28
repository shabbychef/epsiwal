# /usr/bin/r
#
# Copyright 2026-2026 Steven E. Pav. All Rights Reserved.
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
# Created: 2026.05.27
# Copyright: Steven E. Pav, 2026
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) { set.seed(as.integer(charToRaw(str))) }
#UNFOLD

context("mle_connorm")
test_that("runs without error",{#FOLDUP
	set.char.seed("fc742d8c-3454-440a-afd6-0bfd9d27e4c5")

	set.seed(1234)
	n <- 10
	y <- rnorm(n)
	A <- matrix(rnorm(n*(n-3)),ncol=n)
	# make sure y satisfies Ay <= b
	b <- A%*%y + runif(nrow(A))
	Sigma <- diag(runif(n))
	mu <- rnorm(n)
	eta <- rnorm(n)

	expect_error(mval <- mle_connorm(y=y,A=A,b=b,eta=eta,Sigma=Sigma),NA)
	expect_true(is.numeric(mval))
	expect_length(mval, 1)
})#UNFOLD

test_that("mle is consistent with ci",{#FOLDUP
	set.char.seed("3e742d8c-3454-440a-afd6-0bfd9d27e4c5")

	set.seed(4321)
	n <- 5
	y <- rnorm(n)
	A <- matrix(rnorm(n*(n-2)),ncol=n)
	b <- A%*%y + runif(nrow(A))
	Sigma <- diag(runif(n))
	eta <- rnorm(n)

	mval <- mle_connorm(y=y,A=A,b=b,eta=eta,Sigma=Sigma)
	
	# The MLE should correspond to some p-value in ci_connorm.
	# Actually, for symmetric distributions it's 0.5.
	# For truncated normal it's not 0.5.
	
	# Check that it satisfies the MLE equation: E_mu[X] = etay
	stp <- psetup(y=y,A=A,b=b,eta=eta,Sigma_eta=Sigma %*% eta)
	
	expected_val <- function(etamu) {
	  sigma <- sqrt(stp$etaSeta)
	  alpha <- (stp$Vminus - etamu) / sigma
	  beta <- (stp$Vplus - etamu) / sigma
	  
	  log_den <- .log_pnorm_diff(alpha, beta)
	  l_phi_alpha <- dnorm(alpha, log=TRUE)
	  l_phi_beta <- dnorm(beta, log=TRUE)
	  
	  if (l_phi_alpha > l_phi_beta) {
	    log_num <- l_phi_alpha + log1p(-exp(l_phi_beta - l_phi_alpha))
	    millsrat <- exp(log_num - log_den)
	  } else if (l_phi_alpha < l_phi_beta) {
	    log_num <- l_phi_beta + log1p(-exp(l_phi_alpha - l_phi_beta))
	    millsrat <- -exp(log_num - log_den)
	  } else {
	    millsrat <- 0
	  }
	  etamu + sigma * millsrat
	}
	
	expect_equal(expected_val(mval), stp$etay, tolerance=1e-6)
})#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
