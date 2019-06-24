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
