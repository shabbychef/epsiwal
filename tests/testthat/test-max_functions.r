
context("max functions")

test_that("mle_connorm_max matches mle_connorm", {
  set.seed(1234)
  yk <- 3
  yk1 <- c(2, 1)
  sigma <- 1
  rho <- 0.5
  
  mle_new <- mle_connorm_max(yk=yk, yk1=yk1, sigma=sigma, rho=rho)
  
  n <- length(yk1) + 1
  y <- c(yk, yk1)
  eta <- c(1, rep(0, n-1))
  Sigma <- matrix(rho * sigma^2, n, n)
  diag(Sigma) <- sigma^2
  A <- matrix(0, n-1, n)
  A[, 1] <- -1
  for (i in 1:(n-1)) A[i, i+1] <- 1
  b <- rep(0, n-1)
  
  mle_old <- mle_connorm(y=y, A=A, b=b, eta=eta, Sigma=Sigma)
  
  expect_equal(mle_new, mle_old, tolerance=1e-7)
})

test_that("ci_connorm_max matches ci_connorm", {
  set.seed(1234)
  yk <- 3
  yk1 <- c(2, 1)
  sigma <- 1
  rho <- 0.5
  
  ci_new <- ci_connorm_max(yk=yk, yk1=yk1, sigma=sigma, rho=rho, level=0.05)
  
  n <- length(yk1) + 1
  y <- c(yk, yk1)
  eta <- c(1, rep(0, n-1))
  Sigma <- matrix(rho * sigma^2, n, n)
  diag(Sigma) <- sigma^2
  A <- matrix(0, n-1, n)
  A[, 1] <- -1
  for (i in 1:(n-1)) A[i, i+1] <- 1
  b <- rep(0, n-1)
  
  ci_old <- ci_connorm(y=y, A=A, b=b, eta=eta, Sigma=Sigma, level=0.05)
  
  expect_equal(ci_new, ci_old, tolerance=1e-7)
})

test_that("mle_connorm_max handles vector yk1", {
  set.seed(4321)
  yk <- 5
  yk1 <- rnorm(5, mean=2, sd=1)
  sigma <- 1.2
  rho <- 0.2
  
  mle_new <- mle_connorm_max(yk=yk, yk1=yk1, sigma=sigma, rho=rho)
  
  n <- length(yk1) + 1
  y <- c(yk, yk1)
  eta <- c(1, rep(0, n-1))
  Sigma <- matrix(rho * sigma^2, n, n)
  diag(Sigma) <- sigma^2
  A <- matrix(0, n-1, n)
  A[, 1] <- -1
  for (i in 1:(n-1)) A[i, i+1] <- 1
  b <- rep(0, n-1)
  
  mle_old <- mle_connorm(y=y, A=A, b=b, eta=eta, Sigma=Sigma)
  
  expect_equal(mle_new, mle_old, tolerance=1e-7)
})

test_that("mle_connorm_max runs without error", {
  set.seed(1234)
  mus <- seq(-1,1,length.out=100)
  ys  <- rnorm(length(mus),mean=mus,sd=1)
  yss <- sort.int(ys, decreasing=TRUE, index.return=TRUE)
  k <- 5
  yk <- yss$x[k]
  yk1 <- yss$x[k+1]
  expect_error(est_muk <- mle_connorm_max(yk, yk1, sigma=1),NA)
})

test_that("mle_connorm_max sensible values", {
  yk1 <- 4
  sigma <- 1.5
  for (mu in c(-1:5)) {
    # E[X] = mu + sigma * phi(alpha) / (1 - Phi(alpha)) where alpha = (yk1 - mu) / sigma
    # This is for one-sided truncation X >= yk1
    yk <- mu + sigma * dnorm((yk1 - mu)/sigma) / pnorm((yk1 - mu)/sigma, lower.tail=FALSE)
    if (yk > yk1) {
      expect_error(est_muk <- mle_connorm_max(yk, yk1, sigma=sigma),NA)
      expect_equal(est_muk, mu, tolerance=1e-6)
    }
  }
})
