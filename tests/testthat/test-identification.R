test_that("recursive covariance and multi-horizon signs survive normalization", {
  s <- jcer_simulate_gvar(n_obs = 60)
  f <- jcer_bgvar(s$data, s$W, draws = 8, seed = 13)
  labels <- f$prepared$labels; w <- f$worlds[[1]]
  ord <- rev(labels)
  P <- matrix(0, length(labels), length(labels))
  for (j in seq_along(labels)) {
    z <- jcer_irf(f, labels[j], labels, horizon = 0, identification = "recursive", ordering = ord)
    P[, j] <- z$draws[, 1, 1]
  }
  expect_equal(tcrossprod(P), w$Sigma_u, tolerance = 1e-12)
  r <- data.frame(variable = c("US.r_short", "US.y", "US.y"), horizon = c(0, 0, 1), sign = c(1, -1, -1))
  z <- jcer_irf_sign(f, r, "US.r_short", c("US.r_short", "US.y"), horizon = 3, seed = 11, max_tries = 3000)
  expect_true(all(z$draws[2, 1:2, ] < 0))
  expect_equal(as.numeric(z$draws[1, 1, ]), rep(.25, 8))
  zz <- jcer_irf_sign(f, r, "US.r_short", c("US.r_short", "US.y"), horizon = 3, seed = 11, max_tries = 3000,
                       normalization_tol = 1e-6)
  expect_equal(z$draws, zz$draws)
})

test_that("disabling current foreign regressors produces identity G0", {
  s <- jcer_simulate_gvar(n_obs = 40)
  f <- jcer_bgvar(s$data, s$W, draws = 4, contemporaneous = FALSE, seed = 2)
  expect_equal(f$worlds[[1]]$G0, diag(6))
  expect_equal(ncol(f$prepared$blocks[[1]]$Z), 5L)
})

test_that("synthetic percent and decimal DGPs represent the same scaled process", {
  a <- jcer_simulate_gvar(n_obs = 30, seed = 3)
  b <- jcer_simulate_gvar(n_obs = 30, seed = 3, rate_unit = "decimal")
  for (i in seq_along(a$data)) {
    expect_equal(a$data[[i]][, "y"], b$data[[i]][, "y"], tolerance = 1e-12)
    expect_equal(a$data[[i]][, "r_short"] / 100, b$data[[i]][, "r_short"], tolerance = 1e-12)
  }
})
