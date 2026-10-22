test_that("installed BGVAR historical-data IRFs match selected-vector processing", {
  skip_if_not_installed("BGVAR")
  if (utils::packageVersion("BGVAR") < "2.7.0") skip("compatibility baseline is BGVAR 2.7.0 or later")
  e <- new.env(parent = emptyenv())
  utils::data("eerData", package = "BGVAR", envir = e)
  cc <- c("US", "JP", "EA", "KR"); vv <- c("y", "Dp", "stir")
  d <- lapply(e$eerData[cc], function(x) {
    a <- x[, vv, drop = FALSE]
    matrix(as.numeric(a), nrow(a), ncol(a), dimnames = list(NULL, vv))
  })
  W <- e$W.trade0012[cc, cc, drop = FALSE]; W <- W / rowSums(W)
  set.seed(412)
  invisible(capture.output(fit <- BGVAR::bgvar(d, W, draws = 10, burnin = 10,
    prior = "MN", SV = FALSE, eigen = FALSE, expert = list(applyfun = lapply), verbose = FALSE)))
  sh <- BGVAR::get_shockinfo("girf"); sh$shock <- "US.stir"; sh$scale <- .0025
  invisible(capture.output(ref <- BGVAR::irf(fit, n.ahead = 3, shockinfo = sh,
    quantiles = c(.05, .16, .5, .84, .95), expert = list(save.store = TRUE), verbose = FALSE)))
  fast <- jcer_irf_bgvar(fit, "US.stir", horizon = 3, shock_size = .0025, rate_unit = "decimal")
  raw <- array(ref$IRF_store[, 1, , ], dim(fast$draws))
  qref <- aperm(array(ref$posterior[, , 1, ], c(dim(raw)[1:2], 5)), c(3, 1, 2))
  expect_lt(max(abs(raw - fast$draws)), 1e-10)
  expect_lt(max(abs(qref - fast$quantiles)), 1e-10)
  expect_equal(as.numeric(fast$draws["US.stir", "0", ]), rep(.0025, 10))
})
