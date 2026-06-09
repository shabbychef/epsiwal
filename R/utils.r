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

# Lee eqns (5.4), (5.5), (5.6)
Vfuncs <- function(z,A,b,ccc) {
  Az <- A %*% z
  Ac <- A %*% ccc
  bres <- b - Az
  brat <- bres
  suppressWarnings({
    brat[Ac!=0] <- brat[Ac!=0] / Ac[Ac!=0]
    Vminus <- max(brat[Ac < 0])
    Vplus  <- min(brat[Ac > 0])
    Vzero  <- min(bres[Ac == 0])
  })
  list(Vminus=Vminus,Vplus=Vplus,Vzero=Vzero)
}
psetup <- function(y,A,b,eta,Sigma_eta) { 
  etaSeta <- as.numeric(t(eta) %*% Sigma_eta)
  ccc <- Sigma_eta / etaSeta  # eqn (5.3)
  etay <- as.numeric(t(eta) %*% y)
  zzz <- y - ccc * etay  # eqn (5.2)
  Vfs <- Vfuncs(zzz,A,b,ccc)
  c(Vfs,list(etay=etay,etaSeta=etaSeta))
}
# assumes we are testing yk is the maximum against yk1
# and Sigma is compound symmetric with variance sigma^2
# and common equicorrelation of rho.
psetup_max <- function(yk,yk1,sigma,rho) {
	list(Vminus=(max(yk1) - rho * yk) / (1 - rho),Vplus=Inf,Vzero=Inf,etay=yk,etaSeta=sigma**2)
}

# returns log(pnorm(b) - pnorm(a))
# assumes a <= b
.log_pnorm_diff <- function(a, b) {
  out <- rep(NA_real_, length(a))
  
  # Case 1: a == b
  idx_eq <- (a == b)
  out[idx_eq] <- -Inf
  
  # Case 2: a + b > 0
  # This includes the case where both are positive (a > 0)
  # and the across-zero case where b is further from zero.
  # We use the survival function Q = 1 - Phi.
  # Phi(b) - Phi(a) = Q(a) - Q(b)
  idx_upper <- !idx_eq & (a + b > 0)
  if (any(idx_upper)) {
    lq1 <- pnorm(a[idx_upper], lower.tail=FALSE, log.p=TRUE)
    lq2 <- pnorm(b[idx_upper], lower.tail=FALSE, log.p=TRUE)
    # logspace_sub(lq1, lq2)
    out[idx_upper] <- lq1 + log1p(-exp(lq2 - lq1))
  }
  
  # Case 3: a + b <= 0
  # This includes the case where both are negative (b < 0)
  # and the across-zero case where a is further from zero.
  # We use the CDF Phi.
  # Phi(b) - Phi(a)
  idx_lower <- !idx_eq & !idx_upper
  if (any(idx_lower)) {
    lp1 <- pnorm(a[idx_lower], lower.tail=TRUE, log.p=TRUE)
    lp2 <- pnorm(b[idx_lower], lower.tail=TRUE, log.p=TRUE)
    # logspace_sub(lp2, lp1)
    out[idx_lower] <- lp2 + log1p(-exp(lp1 - lp2))
  }
  
  out
}

# returns log(exp(a) + exp(b))
.logspace_add <- function(a, b) {
  pmax(a, b) + log1p(exp(-abs(a - b)))
}

# returns the Mill's ratio Q(x)/phi(x) stably using a continued fraction.
# Note that this is only accurate for large x (e.g. x > 5).
# For small x, more terms are needed, or one should use
# exp(pnorm(x,lower.tail=FALSE,log.p=TRUE) - dnorm(x,log=TRUE))
.mills_ratio <- function(x, n=20) {
  if (is.infinite(x)) return(0)
  res <- 0
  for (i in n:1) {
    res <- i / (x + res)
  }
  1 / (x + res)
}

# returns the "generalized inverse Mill's ratio" (phi(a) - phi(b)) / (Phi(b) - Phi(a))
# which is the adjustment term for the mean of a doubly-truncated normal.
.inv_mills_diff <- function(a, b) {
  # if both are large positive, we use the Mill's ratio to avoid subtraction of small values
  if (a > 10 && b > a) {
    L <- -0.5 * (b^2 - a^2)
    # inv_mills_diff = (phi(a) - phi(b)) / (Q(a) - Q(b)) 
    #                = (1 - exp(L)) / (r(a) - exp(L) r(b))
    eL <- exp(L)
    ra <- .mills_ratio(a)
    rb <- .mills_ratio(b)
    num <- - expm1(L)
    den <- ra - eL * rb
    return(num / den)
  }
  # if both are large negative, exploit symmetry
  if (b < -10 && a < b) {
    return(-.inv_mills_diff(-b, -a))
  }
  
  # default case using logs
  log_den <- .log_pnorm_diff(a, b)
  l_phi_a <- dnorm(a, log=TRUE)
  l_phi_b <- dnorm(b, log=TRUE)
  
  if (l_phi_a > l_phi_b) {
    log_num <- l_phi_a + log(-expm1(l_phi_b - l_phi_a))
    return(exp(log_num - log_den))
  } else if (l_phi_a < l_phi_b) {
    log_num <- l_phi_b + log(-expm1(l_phi_a - l_phi_b))
    return(-exp(log_num - log_den))
  } 
  return(0)
}

#' @importFrom stats uniroot
## invert the ptn function to find y at a given pval.
#qtn <- function(p,A,b,eta,mu,Sigma=NULL,
                #Sigma_eta=Sigma %*% eta,eta_mu=as.numeric(t(eta) %*% mu),
                #intvl=c(-10,10),lower.tail=TRUE) {
  #if (! lower.tail) { p <- 1 - p }

  #f <- function(y) {
    ## rather than this, which is numerically instable
    ##ptn(y,A,b,eta=eta,mu=mu,Sigma=Sigma,etamu=etamu) - p
    ## this
    #stp <- psetup(y=y,A=A,b=b,eta=eta,Sigma_eta=Sigma_eta)
    #phis <- pnorm(c(stp$etay,stp$Vminus,stp$Vplus),mean=eta_mu,sd=sqrt(stp$etaSeta))
    ##(phis[1] - phis[2]) - p * (phis[3] - phis[2])
    #phis[1] + (p-1) * phis[2] - p * phis[3]
  #}
  ## oh, yeah this is a hack
  #trypnts <- seq(from=min(intvl),to=max(intvl),length.out=101)
  #ys <- sapply(trypnts,f)
  #dsy <- diff(sign(ys))
  #if (any(dsy < 0)) {
    #widx <- which(dsy < 0)
    #intvl <- trypnts[widx + c(0,1)]
  #}
  #uniroot(f=f,interval=intvl,extendInt='yes')$root
#}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
