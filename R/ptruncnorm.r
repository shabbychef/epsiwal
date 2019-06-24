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
# Created: 2019.06.12
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @title ptruncnorm .
#'
#' @description 
#'
#' Cumulative distribution of the truncated normal function.
#'
#' @param a  vector of the left truncation value(s).
#' @param b  vector of the right truncation value(s).
#' @inheritParams stats::pnorm
#' @return The distribution function of the truncated normal.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @note Input are recycled as possible.
#' @template etc
#' @importFrom stats pnorm
#' @references
#'
#' Hattaway, James T. "Parameter estimation and hypothesis testing for the truncated normal distribution 
#' with applications to introductory statistics grades." BYU Masters Thesis (2010).
#' \url{https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=3052&context=etd}
#'
#' @examples 
#' y <- ptruncnorm(seq(-5,5,length.out=101), a=-1, b=2)
#' @export 
ptruncnorm <- function(q, mean=0, sd=1, a=-Inf, b=Inf, lower.tail=TRUE, log.p=FALSE) {
	phiq <- pnorm(pmin(pmax(q,a),b),mean=mean,sd=sd)
	phia <- pnorm(a,mean=mean,sd=sd)
	phib <- pnorm(b,mean=mean,sd=sd)
	deno <- phib - phia

	if (lower.tail) {
		ret <- (phiq - phia) / deno
	} else {
		ret <- (phib - phiq) / deno
	}

	if (log.p) { ret <- log(ret) }
	ret
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
