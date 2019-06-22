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
# Author: Steven E. Pav <steven@gilgamath.com>
# Comments: Steven E. Pav

# helpers#FOLDUP
set.char.seed <- function(str) { set.seed(as.integer(charToRaw(str))) }
#UNFOLD

context("ptruncnorm")
test_that("runs without error",{#FOLDUP
	set.char.seed("d106bf01-7ad5-49b2-8715-bb8148a0e294")

	expect_error(ptruncnorm(1,a=0,b=2),NA)
})#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
