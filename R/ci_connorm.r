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
# Created: 2019.06.22
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @title ci_connorm .
#'
#' @description 
#'
#' Confidence intervals on normal mean, subject to linear constraints.
#'
#' @details
#'
#' Inverts the constrained normal inference procedure described
#' by Lee \emph{et al.}
#'
#' Let \eqn{y} be multivariate normal with unknown mean \eqn{\mu}
#' and known covariance \eqn{\Sigma}. Conditional on \eqn{Ay \le b}{Ay <= b}
#' for conformable matrix \eqn{A} and vector \eqn{b}, and given
#' constrast vector \eqn{eta} and level \eqn{p}, we compute
#' \eqn{\eta^{\top}\mu} such that the cumulative distribution of
#' \eqn{\eta^{\top}y} equals \eqn{p}.
#'
#' @param p  a vector of probabilities for which we return
#' equivalent \eqn{\eta^{\top}\mu}.
#' @param level  if \code{p} is not given, we set it by default to
#' \code{c(level/2,1-level/2)}.
#' @inheritParams pconnorm
#' @return The values of \eqn{\eta^{\top}\mu} which have the corresponding
#' CDF.
#' @note 
#' An error will be thrown if we do not observe \eqn{A y \le b}{A y <= b}.
#' @seealso the CDF function, \code{pconnorm}.
#' @template etc
#' @template ref-lee
#' @importFrom stats uniroot
#' @examples
#' set.seed(1234)
#' n <- 10
#' y <- rnorm(n)
#' A <- matrix(rnorm(n*(n-3)),ncol=n)
#' b <- A%*%y + runif(nrow(A))
#' Sigma <- diag(runif(n))
#' mu <- rnorm(n)
#' eta <- rnorm(n)
#' 
#' pval <- pconnorm(y=y,A=A,b=b,eta=eta,mu=mu,Sigma=Sigma)
#' cival <- ci_connorm(y=y,A=A,b=b,eta=eta,Sigma=Sigma,p=pval)
#' stopifnot(abs(cival - sum(eta*mu)) < 1e-4)
#'
#' @export
ci_connorm <- function(y,A,b,eta,Sigma=NULL,p=c(level/2,1-(level/2)),
											 level=0.05,Sigma_eta=Sigma %*% eta) {

	stp <- psetup(y=y,A=A,b=b,eta=eta,Sigma_eta=Sigma_eta)
	stp$sigma <- sqrt(stp$etaSeta)
	# as a hack, a sane range of eta'mu is eta'y +/- 5 sigma
	rang <- stp$etay + 5 * c(-1,1) * stp$sigma

	# you want this, but there are numerical issues: 
	#f <- function(etamu,ap) { F_fnc(x=etay,a=Vfs$Vminus,b=Vfs$Vplus,mu=etamu,sigmasq=etaSeta) - ap } 
	f <- function(etamu,ap) { 
		phis <- pnorm(c(stp$etay,stp$Vminus,stp$Vplus),mean=etamu,sd=stp$sigma)
		#(phis[1] - phis[2]) - p * (phis[3] - phis[2])
		phis[1] + (ap-1) * phis[2] - ap * phis[3]
	}

	sp <- sort.int(p,index.return=TRUE)
	resu <- rep(NA,length(sp$x))

	for (lll in (1:length(sp$x))) {
		nextp <- sp$x[lll]
		if (nextp==0) {
			rootval <- Inf
		} else if (nextp==1) {
			rootval <- -Inf
		} else {
			trypnts <- seq(rang[1],rang[2],length.out=101)
			ys <- sapply(trypnts,f,ap=nextp)
			dsy <- diff(sign(ys))
			intvl <- rang
			if (any(dsy < 0)) {
				widx <- which(dsy < 0)
				intvl <- trypnts[widx + c(0,1)]
			} else {
				delr <- rang[2] - rang[1]
				rang[1] <- rang[1] - 2 * delr
				trypnts <- seq(rang[1],rang[2],length.out=101)
				ys <- sapply(trypnts,f,ap=nextp)
				dsy <- diff(sign(ys))
				if (any(dsy < 0)) {
					widx <- which(dsy < 0)
					intvl <- trypnts[widx + c(0,1)]
				}
			}
			rootval <- uniroot(f=f,interval=intvl,extendInt='yes',ap=nextp)$root
		}
		resu[sp$ix[lll]] <- rootval
		# fix rang
		if (!is.infinite(rootval)) { 
			rang[1] <- rootval
		}
	}
	resu
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
