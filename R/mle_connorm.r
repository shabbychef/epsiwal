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

#' @title mle_connorm .
#'
#' @description 
#'
#' Maximum likelihood estimate of normal mean, subject to linear constraints.
#'
#' @details
#'
#' Computes the maximum likelihood estimate of unknown mean of a normal vector
#' conditional on linear constraints.
#'
#' Let \eqn{y} be multivariate normal with unknown mean \eqn{\mu}
#' and known covariance \eqn{\Sigma}. Conditional on \eqn{Ay \le b}{Ay <= b}
#' for conformable matrix \eqn{A} and vector \eqn{b}, and given
#' constrast vector \eqn{eta}, we compute
#' the maximum likelihood estimate of \eqn{\eta^{\top}\mu}.
#'
#' @param ... dots are passed to \code{uniroot}.
#' @inheritParams pconnorm
#' @return The maximum likelihood estimate of \eqn{\eta^{\top}\mu}.
#' @seealso the confidence interval function, \code{\link{ci_connorm}},
#' the CDF function, \code{\link{pconnorm}},
#' the special case code for conditioning on the max, \code{\link{mle_connorm_max}}
#' @template etc
#' @template ref-reid
#' @importFrom stats uniroot dnorm pnorm
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
#' mval <- mle_connorm(y=y,A=A,b=b,eta=eta,Sigma=Sigma)
#' # try again, but control tolerance:
#' mval <- mle_connorm(y=y,A=A,b=b,eta=eta,Sigma=Sigma,tol=1e-8)
#'
#' @export
mle_connorm <- function(y,A,b,eta,Sigma=NULL,
                        Sigma_eta=Sigma %*% eta,
                        ...) {

  stp <- psetup(y=y,A=A,b=b,eta=eta,Sigma_eta=Sigma_eta)
  .mle_connorm_core(etay=stp$etay, sigma=sqrt(stp$etaSeta), 
										Vminus=stp$Vminus, Vplus=stp$Vplus, ...)
}

#' @title mle_connorm_max .
#'
#' @description 
#'
#' Maximum likelihood estimate of normal mean, conditioning on the max.
#'
#' @details
#'
#' Computes the maximum likelihood estimate of unknown mean of a normal vector
#' conditional on the one element being the maximum.
#'
#' Let \eqn{y} be multivariate normal with unknown mean \eqn{\mu}
#' and known covariance \eqn{\Sigma}. We assume that \eqn{\Sigma}
#' is compound symmetric with common variance \eqn{\sigma^2} and 
#' common correlation \eqn{\rho}. 
#'
#' Conditional on \eqn{y_k \ge y_i}{y_k >= y_i} for all \eqn{i},
#' we compute the maximum likelihood estimate of \eqn{\mu_k}.
#'
#' @param yk the observed maximum value, \eqn{y_k}.
#' @param yk1 a vector of the other observed values, \eqn{y_{k1}}, or just the
#' scalar second largest value.
#' @param sigma the common standard deviation.
#' @param rho the common correlation.
#' @param ... dots are passed to \code{uniroot}.
#' @return The maximum likelihood estimate of \eqn{\mu_k}.
#' @seealso the confidence interval function, \code{\link{ci_connorm_max}},
#' the CDF function, \code{\link{pconnorm}},
#' the more general version, \code{\link{mle_connorm}}.
#' @template etc
#' @template ref-reid
#' @importFrom stats uniroot
#' @export
mle_connorm_max <- function(yk, yk1, sigma=1.0, rho=0, ...) {
  # (Ac)_i = rho - 1 < 0
  # (Az)_i = yi - rho * yk
  # b = 0
  # theta >= (b - Az) / Ac = (0 - (yi - rho * yk)) / (rho - 1) = (yi - rho * yk) / (1 - rho)
  stp <- psetup_max(yk=yk,yk1=yk1,sigma=sigma,rho=rho)
  
  .mle_connorm_core(etay=yk, sigma=sigma, 
										Vminus=stp$Vminus, Vplus=stp$Vplus, ...)
}

.mle_connorm_core <- function(etay, sigma, Vminus, Vplus, ...) {
  alpha0 <- (Vminus - etay) / sigma
  beta0 <- (Vplus - etay) / sigma

  # the MLE satisfies E[X] = etay
  # Let theta = (etamu - etay) / sigma
  # then etamu = etay + sigma * theta
  # alpha = (Vminus - (etay + sigma * theta)) / sigma = alpha0 - theta
  # beta = (Vplus - (etay + sigma * theta)) / sigma = beta0 - theta
  # E[X] = etay + sigma * (theta + inv_mills_diff(alpha, beta))
  # Solve E[X] = etay  =>  theta + inv_mills_diff(alpha0 - theta, beta0 - theta) = 0
  
  f <- function(theta) {
    alpha <- alpha0 - theta
    beta <- beta0 - theta
    
    # E[X] = etamu + sigma * (phi(alpha) - phi(beta)) / (Phi(beta) - Phi(alpha))
    # Solve E[X] = etay  =>  theta + inv_mills_diff(alpha, beta) = 0
    
    theta + .inv_mills_diff(alpha, beta)
  }

  # E[X] is in (Vminus, Vplus). etay must be in (Vminus, Vplus).
  # theta = (etamu - etay) / sigma.
  # If etay -> Vminus, E[X] -> Vminus, so theta -> -Inf.
  # If etay -> Vplus, E[X] -> Vplus, so theta -> Inf.
  # Bounds on theta: (Vminus - etay)/sigma < theta < (Vplus - etay)/sigma
  
  # Search for interval
  # Use the constraints to bound theta
  lower <- alpha0 
  upper <- beta0
  
  # If etay is extremely close to boundary, clamp
  if (is.finite(lower) && abs(lower) < 1e-9) lower <- -1
  if (is.finite(upper) && abs(upper) < 1e-9) upper <- 1

  # If alpha0 or beta0 are infinite or extreme, clamp them
  lower <- max(-100, lower)
  upper <- min(100, upper)
  
  # Ensure interval is valid
  if (lower >= upper) {
    if (lower > 0) lower <- lower - 1
    if (upper < 0) upper <- upper + 1
  }
  
  dots <- list(...)
  if (!("tol" %in% names(dots))) {
    dots$tol <- .Machine$double.eps^0.75
  }
  
  args <- c(list(f=f, interval=c(lower, upper), extendInt='yes'), dots)
  theta_hat <- do.call(uniroot, args)$root
  etay + sigma * theta_hat
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
