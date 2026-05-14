# /usr/bin/r
#
# Copyright 2019-2026 Steven E. Pav. All Rights Reserved.
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
# Created: 2026.05.13
# Copyright: Steven E. Pav, 2026
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) { set.seed(as.integer(charToRaw(str))) }
#UNFOLD

context("truncnorm_mle")
test_that("runs without error",{#FOLDUP
	set.char.seed("e24ca934-cb9c-421b-a424-a9dc63d4b077")

	set.seed(1234)
	mus <- seq(-1,1,length.out=100)
	ys  <- rnorm(length(mus),mean=mus,sd=1)
	yss <- sort.int(ys, decreasing=TRUE, index.return=TRUE)
	k <- 5
	yk <- yss$x[k]
	yk1 <- yss$x[k+1]
	expect_error(est_muk <- truncnorm_mle(yk, yk1, sigma=1),NA)
})#UNFOLD
test_that("sensible values",{#FOLDUP
	set.char.seed("29bd6f7f-4bf5-4d2e-a587-0c8f689a0f0b")
	yk1 <- 4
	sigma <- 1.5
	for (mu in c(-1:5)) {
		yk <- mu + sigma * (dnorm((abs(yk1) - mu)/sigma) - dnorm((-abs(yk1) - mu)/sigma)) / (1 + pnorm((-abs(yk1) - mu)/sigma) - pnorm((abs(yk1) - mu)/sigma))
		if (abs(yk) > abs(yk1)) {
			expect_error(est_muk <- truncnorm_mle(yk, yk1, sigma=sigma),NA)
			# need better tolerances from uniroot?
			expect_equal(est_muk, mu, tolerance=1e-6)
		}
	}
})#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
