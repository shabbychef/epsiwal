context("Numerical stability in tails")

test_that("ci_connorm handles extreme tails", {
  # Simulation setup from bug report
  set.seed(123)
  p <- 100
  srs <- rnorm(p)
  sigma <- 1.2
  mySigma <- diag(sigma, p)

  # Define polyhedral constraints for the maximum element
  yord <- order(srs, decreasing = TRUE)
  revo <- seq_len(p)
  revo[yord] <- revo
  A <- cbind(-1, diag(p - 1))[, revo]
  nu <- rep(0, p)
  nu[yord[1]] <- 1
  b <- rep(0, p - 1)

  # Compute polyhedral median
  # This used to fail by returning a value around -7.12
  val_epsiwal <- epsiwal::ci_connorm(y = srs, A = A, b = b, eta = nu, Sigma = mySigma, p = 0.5)

  # Verify the median property: P(eta'y <= eta'y_obs | Ay <= b) should be 0.5
  x <- srs[yord[1]]
  a <- srs[yord[2]]
  sd_eta <- sqrt(sigma)

  # Stable probability calculation using log-survival functions
  lQx <- pnorm(x, val_epsiwal, sd_eta, lower.tail = FALSE, log.p = TRUE)
  lQa <- pnorm(a, val_epsiwal, sd_eta, lower.tail = FALSE, log.p = TRUE)
  prob_actual <- 1 - exp(lQx - lQa)

  expect_equal(prob_actual, 0.5, tolerance = 1e-5)
  expect_equal(val_epsiwal, -43.057119, tolerance = 1e-5)
})
