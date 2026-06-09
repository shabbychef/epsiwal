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
#' @seealso the CDF function, \code{\link{pconnorm}}, the MLE function, \code{\link{mle_connorm}},
#' the special case code for conditioning on the max, \code{\link{ci_connorm_max}}
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
  .ci_connorm_core(etay=stp$etay, sigma=sqrt(stp$etaSeta), 
                   Vminus=stp$Vminus, Vplus=stp$Vplus, 
                   p=p, level=level)
}

#' @title ci_connorm_max .
#'
#' @description 
#'
#' Confidence intervals on normal mean, conditioning on the max.
#'
#' @details
#'
#' Computes the confidence interval of unknown mean of a normal vector
#' conditional on the one element being the maximum.
#'
#' Let \eqn{y} be multivariate normal with unknown mean \eqn{\mu}
#' and known covariance \eqn{\Sigma}. We assume that \eqn{\Sigma}
#' is compound symmetric with common variance \eqn{\sigma^2} and 
#' common correlation \eqn{\rho}. 
#'
#' Conditional on \eqn{y_k \ge y_i}{y_k >= y_i} for all \eqn{i},
#' we compute the confidence interval of \eqn{\mu_k}.
#'
#' @param yk the observed maximum value, \eqn{y_k}.
#' @param yk1 a vector of the other observed values, \eqn{y_{k1}}, or just the
#' scalar second largest value.
#' @param sigma the common standard deviation.
#' @param rho the common correlation.
#' @inheritParams ci_connorm
#' @return The values of \eqn{\mu_k} which have the corresponding
#' CDF.
#' @seealso the CDF function, \code{\link{pconnorm}}, the MLE function, \code{\link{mle_connorm_max}},
#' the more general version, \code{\link{ci_connorm}}.
#' @template etc
#' @template ref-lee
#' @importFrom stats uniroot
#' @export
ci_connorm_max <- function(yk, yk1, sigma=1.0, rho=0, p=c(level/2,1-(level/2)),
                           level=0.05) {
  stp <- psetup_max(yk=yk,yk1=yk1,sigma=sigma,rho=rho)
  .ci_connorm_core(etay=yk, sigma=sigma, 
                   Vminus=stp$Vminus, Vplus=stp$Vplus, 
                   p=p, level=level)
}

# does the work of the CI functions
.ci_connorm_core <- function(etay, sigma, Vminus, Vplus, p=c(level/2,1-(level/2)),
                             level=0.05) {

  # as a hack, a sane range of eta'mu is eta'y +/- 5 sigma
  rang <- etay + 5 * c(-1,1) * sigma

  # you want this, but there are numerical issues: 
  #f <- function(etamu,ap) { F_fnc(x=etay,a=Vfs$Vminus,b=Vfs$Vplus,mu=etamu,sigmasq=etaSeta) - ap } 
  f <- function(etamu,ap) { 
    if (ap < 0.5) {
      ptruncnorm(q=etay,mean=etamu,sd=sigma,a=Vminus,b=Vplus,log.p=TRUE) - log(ap)
    } else {
      - (ptruncnorm(q=etay,mean=etamu,sd=sigma,a=Vminus,b=Vplus,lower.tail=FALSE,log.p=TRUE) - log1p(-ap))
    }
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
