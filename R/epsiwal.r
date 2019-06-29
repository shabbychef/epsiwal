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

# Created: 2019-06-10
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' Exact Post Selection Inference with Applications to the Lasso.
#' 
#' This simple package supports the simple procedure outlined in 
#' Lee \emph{et al.} where one observes a normal random variable,
#' then performs inference conditional on some linear inequalities.
#'
#' Suppose \eqn{y} is multivariate normal with mean \eqn{\mu}
#' and covariance \eqn{\Sigma}. Conditional on \eqn{Ay \le b}{Ay <= b},
#' one can perform inference on \eqn{\eta^{\top}\mu}{eta'mu} by
#' transforming \eqn{y} to a truncated normal. 
#' Similarly one can invert this procedure and find confidence intervals on
#' \eqn{\eta^{\top}\mu}{eta'mu}.
#'
#' @section Legal Mumbo Jumbo:
#'
#' epsiwal is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU Lesser General Public License for more details.
#'
#' @template etc
#'
#' @template ref-lee
#' @template ref-pav
#'
#' @name epsiwal
#' @rdname epsiwal
#' @docType package
#' @title Exact Post Selection Inference with Applications to the Lasso.
#' @keywords package
#' @note
#' 
#' This package is maintained as a hobby. 
#'
NULL

#' @title News for package 'epsiwal':
#'
#' @description 
#'
#' News for package \sQuote{epsiwal}
#'
#' \newcommand{\CRANpkg}{\href{https://cran.r-project.org/package=#1}{\pkg{#1}}}
#' \newcommand{\epsiwal}{\CRANpkg{epsiwal}}
#'
#' @section \epsiwal{} Initial Version 0.1.0 (2019-06-28) :
#' \itemize{
#' \item first CRAN release.
#' }
#'
#' @name epsiwal-NEWS
#' @rdname NEWS
NULL

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
