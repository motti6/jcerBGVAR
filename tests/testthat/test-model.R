test_that("weights, sample indices and heterogeneous unsupported input are checked", {
  s <- jcer_simulate_gvar(n_obs = 30)
  expect_equal(jcer_prepare_data(s$data, s$W[c(3, 1, 2), c(2, 3, 1)])$X, do.call(cbind, unname(s$data)), ignore_attr = TRUE)
  bad <- s$W; bad[1, 1] <- .1
  expect_error(jcer_prepare_data(s$data, bad), "self weights")
  bad <- s$W; bad[1, ] <- 0
  expect_error(jcer_prepare_data(s$data, bad), "zero foreign")
  expect_error(jcer_prepare_data(s$data, s$W * 2), "row sums")
  expect_equal(jcer_prepare_data(s$data, s$W * 2, normalize_weights = TRUE)$weights[[1]]$W, s$W)
  d <- s$data; d[[1]][1, 1] <- NA
  expect_error(jcer_prepare_data(d, s$W), "missing")
  d <- s$data; d[[1]] <- d[[1]][, 1, drop = FALSE]
  expect_error(jcer_prepare_data(d, s$W), "identical variable")
  d <- s$data; rownames(d[[1]]) <- as.character(1:30)
  expect_error(jcer_prepare_data(d, s$W), "sample indices")
  expect_error(jcer_prepare_data(s$data, s$W, lags = 27), "too short")
  vw <- setNames(rep(list(s$W), 2), c("y", "r_short"))
  expect_equal(jcer_prepare_data(s$data, vw)$blocks, jcer_prepare_data(s$data, s$W)$blocks)
  # Excluding a country requires explicit normalization and is recorded.
  expect_error(jcer_prepare_data(s$data[1:2], s$W), "row sums")
  p <- jcer_prepare_data(s$data[1:2], s$W, normalize_weights = TRUE)
  expect_identical(p$weights[[1]]$excluded, "EA")
  expect_true(p$weights[[1]]$renormalized)
})

test_that("world maps agree with hand constructed two country model", {
  s <- jcer_simulate_gvar(2, 30, variables = "y")
  p <- jcer_prepare_data(s$data, s$W)
  b <- list(matrix(c(0, .3, .2, .1), 4, 1), matrix(c(0, .4, .1, .05), 4, 1))
  w <- preJCERBGVAR:::.world_system(b, list(matrix(1), matrix(2)), p)
  expect_equal(w$G0, matrix(c(1, -.1, -.2, 1), 2))
  expect_equal(w$G[[1]], matrix(c(.3, .05, .1, .4), 2))
  gi <- solve(w$G0)
  expect_equal(w$Sigma_u, gi %*% diag(c(1, 2)) %*% t(gi))
  b[[1]][3, 1] <- 1; b[[2]][3, 1] <- 1
  expect_error(preJCERBGVAR:::.world_system(b, list(matrix(1), matrix(2)), p), "singular")
  b[[1]][3, 1] <- 0; b[[2]][3, 1] <- 0; b[[1]][2, 1] <- 1.2
  expect_false(preJCERBGVAR:::.world_system(b, list(matrix(1), matrix(2)), p)$stable)
})

test_that("three country two variable coefficients and multi-lag recurrence agree", {
  s <- jcer_simulate_gvar(3, 40, lags = 2)
  p <- jcer_prepare_data(s$data, s$W, lags = 2)
  w <- s$truth$world
  expect_equal(w$G0, diag(6) - kronecker(s$W, diag(.08, 2)))
  A <- matrix(c(.2, 0, -.03, .2), 2)
  expect_equal(w$G[[1]], kronecker(diag(3), A) + kronecker(s$W, diag(.025, 2)))
  impact <- 1:6 / 10
  got <- preJCERBGVAR:::.recur_irf(w, impact, 8)
  ref <- matrix(0, 6, 9); ref[, 1] <- impact
  gi <- solve(w$G0)
  for (h in 1:8) for (l in 1:min(h, 2)) ref[, h + 1] <- ref[, h + 1] + gi %*% w$G[[l]] %*% ref[, h - l + 1]
  expect_equal(got, ref, tolerance = 1e-12)
  expect_true(w$stable)
})

test_that("posterior matches independently evaluated formula and moments", {
  set.seed(17); Z <- cbind(1, matrix(rnorm(90), 30, 3)); Y <- matrix(rnorm(60), 30, 2)
  pr <- preJCERBGVAR:::.prior_mniw(4L, 2L, list())
  po <- preJCERBGVAR:::.mniw_posterior(Y, Z, pr)
  V <- solve(solve(pr$V0) + crossprod(Z))
  B <- V %*% (solve(pr$V0) %*% pr$B0 + crossprod(Z, Y))
  S <- pr$S0 + crossprod(Y) + crossprod(pr$B0, solve(pr$V0) %*% pr$B0) - crossprod(B, solve(V) %*% B)
  expect_equal(po$Bn, B); expect_equal(po$Vn, V); expect_equal(po$Sn, S); expect_equal(po$nun, pr$nu0 + 30)
  dd <- preJCERBGVAR:::.country_draws(list(Y = Y, Z = Z, prior = pr, seed = 9, draws = 3000))
  expect_equal(apply(dd$B, 1:2, mean), B, tolerance = .025)
  expect_equal(apply(dd$Sigma, 1:2, mean), S / (po$nun - 3), tolerance = .035)
  expect_error(preJCERBGVAR:::.prior_mniw(4L, 2L, list(nu0 = 1)), "nu0")
  expect_error(preJCERBGVAR:::.prior_mniw(4L, 2L, list(S0 = diag(c(1, -1)))), "positive definite")
})

test_that("draws are reproducible, prefix invariant, and RNG state is restored", {
  s <- jcer_simulate_gvar(n_obs = 40)
  set.seed(42); old <- .Random.seed
  f <- jcer_bgvar(s$data, s$W, draws = 6, seed = 3)
  expect_identical(old, .Random.seed)
  g <- jcer_bgvar(s$data, s$W, draws = 10, seed = 3)
  expect_identical(f$countries[[1]]$B, g$countries[[1]]$B[, , 1:6])
  expect_identical(jcer_diagnostics(f)$burn_in, 0L)
  expect_error(jcer_bgvar(s$data, s$W, method = "vb"), "validated")
})

test_that("impact units, sign restrictions, Haar covariance and plot labels hold", {
  s <- jcer_simulate_gvar(n_obs = 60)
  f <- jcer_bgvar(s$data, s$W, draws = 15, seed = 3)
  a <- jcer_irf(f, "US.r_short", c("US.r_short", "JP.y"), horizon = 3, shock_bp = 25)
  expect_equal(as.numeric(a$draws[1, 1, ]), rep(.25, 15))
  f$prepared$rate_unit <- "decimal"
  b <- jcer_irf(f, "US.r_short", "US.r_short", horizon = 0, shock_bp = 25)
  expect_equal(as.numeric(b$draws), rep(.0025, 15))
  expect_error(jcer_irf(f, "US.r_short", "JP.y", identification = "recursive"), "ordering")
  q <- preJCERBGVAR:::.haar(6); P <- t(chol(f$worlds[[1]]$Sigma_u))
  expect_equal(tcrossprod(q), diag(6), tolerance = 1e-12)
  expect_equal(tcrossprod(P %*% q), f$worlds[[1]]$Sigma_u, tolerance = 1e-12)
  r <- data.frame(variable = "US.r_short", horizon = 0, sign = 1)
  z <- jcer_irf_sign(f, r, "US.r_short", "US.r_short", horizon = 3, seed = 4)
  expect_equal(as.numeric(z$draws[1, 1, ]), rep(.0025, 15))
  expect_true(all(z$accepted == 1)); expect_true(z$acceptance_rate > 0)
  expect_error(jcer_irf_sign(f, transform(r, sign = -1), "US.r_short", "JP.y", max_tries = 2, seed = 1), class = "JCERBGVAR_sign_failure")
  expect_error(jcer_irf_sign(f, r, "US.r_short", "JP.y", normalization_tol = 1e6, max_tries = 2), "insufficient")
  tf <- tempfile(fileext = ".pdf"); grDevices::pdf(tf)
  expect_invisible(plot(a)); grDevices::dev.off(); unlink(tf)
  expect_true(a$synthetic); expect_match(a$interpretation, "not an identified")
})
