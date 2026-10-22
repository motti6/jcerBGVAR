Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
# Dense validation reference: invert G0 once, propagate all impulse columns,
# then project onto the same selected innovation. No repeated inverse strawman.
dense_reference <- function(w, impact, horizon) {
  gi <- solve(w$G0); A <- lapply(w$G, function(g) gi %*% g)
  m <- nrow(gi); psi <- array(0, c(m, m, horizon + 1L)); psi[, , 1L] <- diag(m)
  for (h in seq_len(horizon)) for (l in seq_len(min(h, length(A))))
    psi[, , h + 1L] <- psi[, , h + 1L] + A[[l]] %*% psi[, , h - l + 1L]
  apply(psi, 3L, function(z) z %*% impact)
}
rows <- list()
for (n in c(3, 12, 20)) {
  s <- jcer_simulate_gvar(n, 120, seed = 1, lags = 2)
  f <- jcer_bgvar(s$data, s$W, lags = 2, draws = 100, seed = 1)
  impacts <- lapply(f$worlds, function(w) w$Sigma_u[, 2] / sqrt(w$Sigma_u[2, 2]))
  for (d in seq_along(f$worlds)) stopifnot(max(abs(dense_reference(f$worlds[[d]], impacts[[d]], 20) -
    preJCERBGVAR:::.recur_irf(f$worlds[[d]], impacts[[d]], 20))) < 1e-11)
  invisible(dense_reference(f$worlds[[1]], impacts[[1]], 20))
  for (r in 1:5) {
    t1 <- system.time(a <- lapply(seq_along(f$worlds), function(d) dense_reference(f$worlds[[d]], impacts[[d]], 20)))[[3]]
    t2 <- system.time(b <- lapply(seq_along(f$worlds), function(d) preJCERBGVAR:::.recur_irf(f$worlds[[d]], impacts[[d]], 20)))[[3]]
    rows[[length(rows) + 1L]] <- data.frame(n_countries = n, n_obs = 120, lags = 2, draws = 100,
      horizon = 20, repetition = r, dense_reference_seconds = t1, selected_vector_seconds = t2,
      speedup = if (t2 > 0) t1 / t2 else NA_real_,
      max_absolute_error = max(abs(unlist(a) - unlist(b))), synthetic = TRUE)
  }
}
ans <- do.call(rbind, rows)
write.csv(ans, "benchmarks/results/irf_algorithm_comparison.csv", row.names = FALSE)
print(aggregate(ans[, c("dense_reference_seconds", "selected_vector_seconds", "speedup", "max_absolute_error")], ans["n_countries"], median))
