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
# Created: 2026.05.13
# Copyright: Steven E. Pav, 2026
# Author: Steven E. Pav <steven@gilgamath.com>
# Comments: Steven E. Pav

#' @title truncnorm_mle computes the MLE truncated norm estimate of population mean.
#'
#' @description 
#'
#' Solves equation (9) of Reid, Taylor and Tibshirani (2014) to estimate the
#' population mean conditional on an estimation procedure.
#'
#' @details
#'
#' Suppose that \eqn{y_i} are independently normally distributed around $\mu_i$ with standard
#' deviation \eqn{\sigma}. Suppose that the \eqn{y_i} are ordered such that
#' \eqn{|y_1| \ge |y_2| \ge \ldots \ge |y_k| \ge |y_{k+1}| \ldots}. This function
#' seeks to estimate \eqn{\mu_k} given \eqn{y_k} and \eqn{y_{k+1}}.
#'
#' @usage
#'
#' truncnorm_mle(yk, yk1, sigma)
#'
#' @param yk  the value of \eqn{y_k}.
#' @param yk1  the value of \eqn{y_{k+1}}. We must have \eqn{|y_k| > |y_{k+1}|}.
#' @param sigma  the (common) standard deviation of noise.
#' @return The estimate of \eqn{\mu_k}.
#' @note 
#' An error will be thrown if we do not observe \eqn{|y_{k}| \ge |y_{k+1}|}.
#' @examples
#'
#' mus <- seq(-1,1,length.out=100)
#' ys  <- rnorm(length(mus),mean=mus,sd=1)
#' yss <- sort.int(ys, decreasing=TRUE, index.return=TRUE)
#' k <- 5
#' yk <- yss$x[k]
#' yk1 <- yss$x[k+1]
#' est_muk <- truncnorm_mle(yk, yk1, sigma=1)
#' # compare to  mus[yss$ix[k]]
#'
#' @template etc
#' @template ref-reid
#' @export
truncnorm_mle <- function(yk, yk1, sigma) {
  if (abs(yk) < abs(yk1)) {
		stop(paste0("expected |yk| >= |yk1|, instead got |",yk,"| < |",yk1,"|."))
	}
	# solvfe this for mu
	fsolve <- function(mu) {
		yplus <- (abs(yk1) - mu) / sigma
		ymins <- (-abs(yk1) - mu) / sigma
	  millsrat <- diff(dnorm(c(ymins, yplus))) / (1 + diff(pnorm(c(yplus, ymins))))
		yk - mu - sigma * millsrat
	}
	# 2FIX: need a better starting estimate.
	interval <- yk + 10*sigma * c(-1,1)
	solv <- uniroot(fsolve, interval, tol=.Machine$double.eps^0.4)
	solv$root
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
