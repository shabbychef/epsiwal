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

#' @title pconnorm .
#'
#' @description 
#'
#' CDF of the conditional normal variate.
#'
#' @details
#'
#' Computes the CDF of the truncated normal conditional
#' on linear constraints, as described in section 5 of Lee \emph{et al.}
#'
#' Let \eqn{y} be multivariate normal with mean \eqn{\mu}
#' and covariance \eqn{\Sigma}. Conditional on \eqn{Ay \le b}{Ay <= b}
#' for conformable matrix \eqn{A} and vector \eqn{b} we compute the
#' CDF of a truncated normal maximally aligned with \eqn{\eta}. 
#' Inference depends on the population parameters only via 
#' \eqn{\eta^{\top}\mu}{eta'mu} and \eqn{\Sigma \eta}{Sigma eta},
#' and only these need to be given.
#'
#' The test statistic is aligned with \eqn{y}, meaning that an output
#' p-value near one casts doubt on the null hypothesis that 
#' \eqn{\eta^{\top}\mu}{eta'mu} is less than the posited value.
#'
#' @param y  an \eqn{n} vector, assumed multivariate normal with mean \eqn{\mu}
#' and covariance \eqn{\Sigma}.
#' @param A  an \eqn{k \times n} matrix of constraints.
#' @param b  a \eqn{k} vector of inequality limits.
#' @param eta  an \eqn{n} vector of the test contrast, \eqn{\eta}.
#' @param mu  an \eqn{n} vector of the population mean, \eqn{\mu}.
#' Not needed if \code{eta_mu} is given.
#' @param Sigma  an \eqn{n \times n} matrix of the population covariance, \eqn{\Sigma}.
#' Not needed if \code{Sigma_eta} is given.
#' @param Sigma_eta  an \eqn{n} vector of \eqn{\Sigma \eta}.
#' @param eta_mu   the scalar \eqn{\eta^{\top}\mu}.
#' @inheritParams stats::pnorm
#' @return The CDF.
#' @note 
#' An error will be thrown if we do not observe \eqn{A y \le b}{A y <= b}.
#' @seealso the confidence interval function, \code{ci_connorm}.
#' @template etc
#' @template ref-lee
#' @export
pconnorm <- function(y,A,b,eta,mu=NULL,Sigma=NULL,
										 Sigma_eta=Sigma %*% eta,eta_mu=as.numeric(t(eta) %*% mu),
										 lower.tail=TRUE,log.p=FALSE) {
	# this is just Lee et. al eqn (5.9)
  stp <- psetup(y=y,A=A,b=b,eta=eta,Sigma_eta=Sigma_eta)
	# ptrunc is Lee et. al eqn (5.8)
  ptruncnorm(q=stp$etay,a=stp$Vminus,b=stp$Vplus,mean=eta_mu,sd=sqrt(stp$etaSeta),
						 lower.tail=lower.tail,log.p=log.p)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
