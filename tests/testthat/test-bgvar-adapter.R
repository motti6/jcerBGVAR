test_that("BGVAR adapter distinguishes country-error and world innovations", {
  gi <- matrix(c(1, .2, .1, 1), 2); se <- diag(c(1, 2))
  F <- matrix(c(.4, .1, -.03, .3), 2)
  labels <- c("US.stir", "JP.y")
  m <- structure(list(xglobal = matrix(0, 10, 2, dimnames = list(NULL, labels)),
    args = list(lags = c(1L, 1L)), stacked.results = list(
      F_large = array(F, c(2, 2, 1, 1)), Ginv_large = array(gi, c(2, 2, 1)),
      S_large = array(se, c(2, 2, 1)), F.eigen = .5)), class = "bgvar")
  a <- jcer_irf_bgvar(m, "US.stir", horizon = 3)
  b <- jcer_irf_bgvar(m, "US.stir", horizon = 3, definition = "world_generalized")
  raw <- gi %*% se[, 1]; impact <- raw * .25 / raw[1]
  expect_equal(a$draws[, 1, 1], as.numeric(impact), ignore_attr = TRUE)
  expect_equal(a$draws[, 4, 1], as.numeric(F %*% F %*% F %*% impact), ignore_attr = TRUE)
  expect_gt(max(abs(a$draws - b$draws)), .01)
  expect_equal(a$draws[1, 1, 1], .25)
  expect_equal(jcer_irf_bgvar(m, "US.stir", horizon = 0, shock_size = .0025, rate_unit = "decimal")$draws[1, 1, 1], .0025)
  expect_error(jcer_irf_bgvar(m, "US.stir", normalization_tol = 2), "near-zero")
  m$stacked.results$F_large[1] <- NA_real_
  expect_error(jcer_irf_bgvar(m, "US.stir"), "nonfinite")
})

test_that("lag-only synthetic DGP and fitting model are aligned", {
  for (p in 1:2) {
    s <- jcer_simulate_gvar(3, 50, lags = p, contemporaneous = FALSE)
    expect_equal(s$truth$world$G0, diag(6))
    expect_match(s$exogeneity, "independent")
    for (l in 1:p) expect_equal(s$truth$world$G[[l]],
      kronecker(diag(3), matrix(c(.4 / p, 0, -.06 / p, .4 / p), 2)) +
        kronecker(s$W, diag(.05 / p, 2)))
    f <- jcer_bgvar(s$data, s$W, draws = 5, lags = p, contemporaneous = FALSE, seed = 1)
    expect_equal(f$worlds[[1]]$G0, diag(6))
  }
  expect_error(jcer_simulate_gvar(persistence = 1.2), "not stable")
})

test_that("BGVAR selected vectors handle multiple lags and singleton selections", {
  A <- matrix(c(.3, .1, -.02, .2), 2); B <- diag(.1, 2)
  F <- array(0, c(2, 2, 2, 2))
  for (d in 1:2) { F[, , 1, d] <- A; F[, , 2, d] <- B }
  labels <- c("US.stir", "JP.y")
  x <- structure(list(xglobal = matrix(0, 10, 2, dimnames = list(NULL, labels)),
    args = list(lags = c(2L, 2L)), stacked.results = list(F_large = F,
      Ginv_large = array(rep(diag(2), 2), c(2, 2, 2)), S_large = array(rep(diag(2), 2), c(2, 2, 2)),
      F.eigen = c(.5, .5))), class = "bgvar")
  got <- jcer_irf_bgvar(x, "US.stir", "JP.y", horizon = 5)
  psi <- array(0, c(2, 2, 6)); psi[, , 1] <- diag(2)
  for (h in 1:5) {
    psi[, , h + 1] <- A %*% psi[, , h]
    if (h > 1) psi[, , h + 1] <- psi[, , h + 1] + B %*% psi[, , h - 1]
  }
  expect_equal(as.numeric(got$draws[1, , 1]), psi[2, 1, ] * .25)
  expect_identical(dim(got$draws), c(1L, 6L, 2L))
})
