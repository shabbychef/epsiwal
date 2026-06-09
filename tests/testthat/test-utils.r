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

context("internal utilities")

test_that(".mills_ratio works", {
  skip_on_cran()
  
  # it is intended for large x, where the continued fraction is stable
  xs <- seq(5, 20, length.out=10)
  for (x in xs) {
    # .mills_ratio(x) should be Q(x) / phi(x)
    expected <- exp(stats::pnorm(x, lower.tail=FALSE, log.p=TRUE) - stats::dnorm(x, log=TRUE))
    expect_equal(.mills_ratio(x), expected, tolerance=1e-7)
  }
  
  # large values
  # for x=10, the naive calculation is still fine
  expected_10 <- exp(stats::pnorm(10, lower.tail=FALSE, log.p=TRUE) - stats::dnorm(10, log=TRUE))
  expect_equal(.mills_ratio(10), expected_10, tolerance=1e-7)
  
  # very large values
  # .mills_ratio(x) approx 1/x for large x
  expect_equal(.mills_ratio(100), 1/100, tolerance=1e-4)
  
  # infinity
  expect_equal(.mills_ratio(Inf), 0)
})

test_that(".logspace_add works", {
  skip_on_cran()
  
  expect_equal(.logspace_add(log(2), log(3)), log(5))
  expect_equal(.logspace_add(100, 100), 100 + log(2))
  expect_equal(.logspace_add(-100, -100), -100 + log(2))
  expect_equal(.logspace_add(100, -100), 100)
})

test_that(".log_pnorm_diff works", {
  skip_on_cran()
  
  # basic case
  expect_equal(.log_pnorm_diff(0, 1), log(pnorm(1) - pnorm(0)))
  
  # symmetry
  expect_equal(.log_pnorm_diff(-1, 0), .log_pnorm_diff(0, 1))
  
  # small diff
  expect_equal(.log_pnorm_diff(10, 10.1), 
               pnorm(10, lower.tail=FALSE, log.p=TRUE) + log1p(-exp(pnorm(10.1, lower.tail=FALSE, log.p=TRUE) - pnorm(10, lower.tail=FALSE, log.p=TRUE))),
               tolerance=1e-12)
  
  # large values
  expect_true(.log_pnorm_diff(10, 11) < 0)
  expect_true(is.finite(.log_pnorm_diff(10, 11)))
  
  # a == b
  expect_equal(.log_pnorm_diff(1, 1), -Inf)
})

test_that(".inv_mills_diff works", {
  skip_on_cran()
  
  # .inv_mills_diff(a, b) should be (phi(a) - phi(b)) / (Phi(b) - Phi(a))
  
  # moderate values
  a <- 0
  b <- 1
  expected <- (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
  expect_equal(.inv_mills_diff(a, b), expected)
  
  # large positive values (uses .mills_ratio path)
  a <- 12
  b <- 13
  # (phi(a) - phi(b)) / (Phi(b) - Phi(a)) = (phi(a) - phi(b)) / (Q(a) - Q(b))
  # divide num and den by phi(a): (1 - phi(b)/phi(a)) / (r(a) - phi(b)/phi(a) * r(b))
  # phi(b)/phi(a) = exp(-0.5 * (b^2 - a^2))
  L <- -0.5 * (b^2 - a^2)
  eL <- exp(L)
  expected_large <- (1 - eL) / (.mills_ratio(a) - eL * .mills_ratio(b))
  expect_equal(.inv_mills_diff(a, b), expected_large)
  
  # large negative values (symmetry path)
  expect_equal(.inv_mills_diff(-13, -12), -.inv_mills_diff(12, 13))
  
  # a == b should return 0 or handles it?
  # current implementation might return NaN or something else if not careful
  # but usually it is called with Vminus < Vplus
})

test_that("Vfuncs works", {
  skip_on_cran()
  
  # Simple case: 1D, A = 1, b = 2, z = 0, ccc = 1
  # Ay <= b => z + ccc*etay <= b => etay <= (b-z)/ccc = 2
  # Ac = 1 (> 0), Vplus = 2
  A <- matrix(1, 1, 1)
  b <- 2
  z <- 0
  ccc <- 1
  vfs <- Vfuncs(z, A, b, ccc)
  expect_equal(vfs$Vplus, 2)
  expect_equal(vfs$Vminus, -Inf)
  expect_equal(vfs$Vzero, Inf)
  
  # Ac < 0 case
  ccc <- -1
  vfs <- Vfuncs(z, A, b, ccc)
  expect_equal(vfs$Vminus, -2)
  expect_equal(vfs$Vplus, Inf)
  
  # Ac = 0 case
  ccc <- 0
  vfs <- Vfuncs(z, A, b, ccc)
  expect_equal(vfs$Vzero, 2)
})

test_that("psetup works", {
  skip_on_cran()
  
  set.seed(123)
  n <- 5
  y <- rnorm(n)
  A <- matrix(rnorm(n*2), ncol=n)
  b <- A %*% y + runif(nrow(A))
  eta <- rnorm(n)
  Sigma_eta <- runif(n)
  
  expect_error(stp <- psetup(y, A, b, eta, Sigma_eta), NA)
  expect_named(stp, c("Vminus", "Vplus", "Vzero", "etay", "etaSeta"))
  
  # check constraint satisfaction
  # Ay <= b => A(z + c*etay) <= b => A*z + A*c*etay <= b
  # if Ac > 0, etay <= (b - Az) / Ac
  # if Ac < 0, etay >= (b - Az) / Ac
  expect_true(stp$etay >= stp$Vminus)
  expect_true(stp$etay <= stp$Vplus)
  expect_true(stp$Vminus <= stp$Vplus)
})

test_that("psetup and psetup_max are equivalent", {
  skip_on_cran()

  set.seed(1234)
  n <- 5
  sigma <- 1.5
  rho <- 0.3
  
  # Equicorrelation matrix
  Sigma <- sigma^2 * ((1 - rho) * diag(n) + rho * matrix(1, n, n))
  
  y <- rnorm(n, mean = 2, sd = sigma) 
  k <- which.max(y)
  yk <- y[k]
  yk1 <- y[-k]
  
  # Setup for psetup
  # A y <= b
  # y_j - y_k <= 0  => (e_j - e_k) y <= 0
  A <- matrix(0, nrow = n - 1, ncol = n)
  j_idx <- (1:n)[-k]
  for (i in seq_along(j_idx)) {
    A[i, j_idx[i]] <- 1
    A[i, k] <- -1
  }
  b <- rep(0, n - 1)
  eta <- rep(0, n)
  eta[k] <- 1
  
  Sigma_eta <- Sigma %*% eta
  
  stp1 <- psetup(y = y, A = A, b = b, eta = eta, Sigma_eta = Sigma_eta)
  stp2 <- psetup_max(yk = yk, yk1 = yk1, sigma = sigma, rho = rho)
  
  expect_equal(stp1$Vminus, stp2$Vminus)
  expect_equal(stp1$Vplus, stp2$Vplus)
  expect_equal(stp1$Vzero, stp2$Vzero)
  expect_equal(stp1$etay, stp2$etay)
  expect_equal(stp1$etaSeta, stp2$etaSeta)
})

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
