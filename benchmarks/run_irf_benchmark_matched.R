Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
outdir <- "benchmarks/results/synthetic_revised_20261022"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Both receive ONLY raw G0/G_l/Sigma_e, and include all setup/factorization,
# impact construction/25bp normalization, propagation and five quantiles.
evaluate <- function(raw, engine, horizon = 20) {
  m <- nrow(raw[[1]]$G0); nd <- length(raw); out <- array(0, c(m, horizon + 1L, nd))
  for (d in seq_len(nd)) {
    w <- raw[[d]]
    if (engine == "dense") {
      gi <- solve(w$G0); F <- lapply(w$G, function(g) gi %*% g)
      su <- gi %*% w$Sigma_e %*% t(gi); impact <- su[, 2] * .25 / su[2, 2]
      psi <- array(0, c(m, m, horizon + 1L)); psi[, , 1] <- diag(m)
      for (h in seq_len(horizon)) for (l in seq_len(min(h, length(F))))
        psi[, , h + 1L] <- psi[, , h + 1L] + F[[l]] %*% psi[, , h - l + 1L]
      out[, , d] <- apply(psi, 3L, function(x) drop(x %*% impact))
    } else {
      fac <- qr(w$G0, LAPACK = TRUE)
      L <- qr.coef(fac, t(chol(w$Sigma_e)))
      col <- drop(L %*% L[2, ]); impact <- col * .25 / col[2]
      ww <- list(G0 = w$G0, G = w$G, factor = fac)
      out[, , d] <- preJCERBGVAR:::.recur_irf(ww, impact, horizon)
    }
  }
  list(draws = out, quantiles = apply(out, c(1, 2), quantile, probs = c(.05, .16, .5, .84, .95), names = FALSE))
}
rows <- list()
for (n in c(3, 12, 20)) {
  s <- jcer_simulate_gvar(n, 120, seed = 1, lags = 2)
  f <- jcer_bgvar(s$data, s$W, lags = 2, draws = 200, seed = 1)
  raw <- lapply(f$worlds, function(w) w[c("G0", "G", "Sigma_e")])
  invisible(evaluate(raw, "dense")); invisible(evaluate(raw, "vector"))
  for (r in 1:5) {
    if (r %% 2) {
      t1 <- system.time(a <- evaluate(raw, "dense"))[[3]]
      t2 <- system.time(b <- evaluate(raw, "vector"))[[3]]
    } else {
      t2 <- system.time(b <- evaluate(raw, "vector"))[[3]]
      t1 <- system.time(a <- evaluate(raw, "dense"))[[3]]
    }
    error <- max(abs(a$draws - b$draws)); stopifnot(error < 1e-10)
    rows[[length(rows) + 1L]] <- data.frame(countries = n, lags = 2, draws = 200, horizon = 20, repetition = r,
      dense_seconds = t1, vector_seconds = t2, speedup = t1 / t2, max_absolute_error = error,
      max_quantile_error = max(abs(a$quantiles - b$quantiles)),
      timing_scope = "raw matrices -> factorization + impact + normalization + all requested responses + five quantiles")
  }
  message("Matched synthetic IRF timing finished: ", n, " countries")
}
ans <- do.call(rbind, rows)
write.csv(ans, file.path(outdir, "matched_irf_timing.csv"), row.names = FALSE)
print(aggregate(ans[, c("dense_seconds", "vector_seconds", "speedup", "max_absolute_error")], ans["countries"], median))
