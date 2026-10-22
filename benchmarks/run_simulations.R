Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
source("benchmarks/reference-gibbs.R")
dir.create("benchmarks/results", recursive = TRUE, showWarnings = FALSE)
full <- "--full" %in% commandArgs(TRUE)
nd <- if (full) 10000 else 2500
s <- jcer_simulate_gvar(seed = 52)
f <- jcer_bgvar(s$data, s$W, draws = nd, seed = 52)
block <- f$prepared$blocks[[1]]; pr <- f$countries[[1]]$posterior$prior
chains <- vector("list", 4); timings <- numeric(4)
for (ch in 1:4) timings[ch] <- system.time(chains[[ch]] <- reference_gibbs(block$Y, block$Z, pr, nd, 1000, 600 + ch))[[3]]
direct <- f$countries[[1]]$B; gibbs <- do.call(cbind, lapply(chains, function(z) matrix(z$B, nrow = length(z$B) / nd)))
dm <- rowMeans(matrix(direct, nrow = length(direct) / nd)); gm <- rowMeans(gibbs)
dv <- apply(matrix(direct, nrow = length(direct) / nd), 1, var); gv <- apply(gibbs, 1, var)
direct_times <- replicate(5, system.time(invisible(preJCERBGVAR:::.country_draws(list(Y = block$Y, Z = block$Z, prior = pr, seed = 52, draws = nd))))[[3]])
di <- median(direct_times)
ess <- sapply(chains, function(z) chain_ess(z$B[2, 1, ]))
rh <- classic_rhat(sapply(chains, function(z) z$B[2, 1, ]))
comparison <- data.frame(comparison = "same conditional MNIW; one country; independent direct vs four-chain Gibbs",
  direct_draws = nd, gibbs_draws_per_chain = nd, burnin_per_chain = 1000,
  direct_seconds = di, gibbs_seconds_per_chain = median(timings),
  speedup_draw_budget = median(timings) / di,
  mean_abs_posterior_mean_difference = mean(abs(dm - gm)), mean_abs_variance_difference = mean(abs(dv - gv)),
  selected_coefficient_Rhat = rh, selected_coefficient_ESS_per_chain = median(ess),
  gibbs_ESS_per_second = median(ess / timings), direct_ESS_per_second = nd / di)
write.csv(comparison, "benchmarks/results/gibbs_comparison.csv", row.names = FALSE)
# Country draws from an independent Gibbs sampler -> full world IRFs, same model.
gfit <- f; gfit$worlds <- vector("list", min(nd, 1000L))
gc <- lapply(seq_along(f$countries), function(i) reference_gibbs(f$prepared$blocks[[i]]$Y,
  f$prepared$blocks[[i]]$Z, f$countries[[i]]$posterior$prior, nd, 1000, 900 + i))
for (d in seq_along(gfit$worlds)) gfit$worlds[[d]] <- preJCERBGVAR:::.world_system(
  lapply(gc, function(x) matrix(x$B[, , d], dim(x$B)[1], dim(x$B)[2])),
  lapply(gc, function(x) matrix(x$Sigma[, , d], dim(x$Sigma)[1], dim(x$Sigma)[2])), f$prepared)
refir <- jcer_irf(gfit, "US.r_short", c("JP.r_short", "JP.y"), horizon = 12, shock_bp = 25)
newir <- jcer_irf(f, "US.r_short", c("JP.r_short", "JP.y"), horizon = 12, shock_bp = 25)
write.csv(jcer_compare_irf(refir, newir), "benchmarks/results/gibbs_irf_comparison.csv", row.names = FALSE)
# Conditional one-step held-out predictive log score: observed foreign paths are
# conditioned on, not jointly forecast. Independent samplers, identical prior.
train <- seq_len(nrow(block$Y) - 10L); holdout <- tail(seq_len(nrow(block$Y)), 10)
dt <- preJCERBGVAR:::.country_draws(list(Y = block$Y[train, , drop = FALSE], Z = block$Z[train, , drop = FALSE],
  prior = pr, seed = 501, draws = nd))
gt <- reference_gibbs(block$Y[train, , drop = FALSE], block$Z[train, , drop = FALSE], pr, nd, 1000, 502)
predictive_score <- function(dr, y, z) {
  logd <- vapply(seq_len(dim(dr$B)[3]), function(d) {
    S <- dr$Sigma[, , d]; U <- chol(S); e <- y - drop(z %*% dr$B[, , d])
    v <- forwardsolve(t(U), e)
    -length(y) / 2 * log(2 * pi) - sum(log(diag(U))) - sum(v^2) / 2
  }, numeric(1))
  ma <- max(logd); ma + log(mean(exp(logd - ma)))
}
scores <- data.frame(holdout_row = holdout,
  direct_log_score = vapply(holdout, function(t) predictive_score(dt, block$Y[t, ], block$Z[t, ]), numeric(1)),
  gibbs_log_score = vapply(holdout, function(t) predictive_score(gt, block$Y[t, ], block$Z[t, ]), numeric(1)),
  target = "conditional one-step density given observed domestic/foreign regressors")
write.csv(scores, "benchmarks/results/predictive_log_scores.csv", row.names = FALSE)
nr <- if (full) 200 else 20
results <- list()
for (rep in seq_len(nr)) {
  sim <- jcer_simulate_gvar(seed = 1000 + rep)
  fit <- jcer_bgvar(sim$data, sim$W, draws = if (full) 1000 else 250, seed = 2000 + rep)
  ir <- jcer_irf(fit, "US.r_short", c("JP.r_short", "JP.y"), horizon = 8)
  truth <- sim$true_irf[match(ir$response, rownames(sim$true_irf)), 1:9, drop = FALSE]
  err <- matrix(ir$quantiles[3, , ], 2) - truth
  results[[rep]] <- data.frame(replication = rep, bias = mean(err), rmse = sqrt(mean(err^2)),
    coverage68 = mean(truth >= ir$quantiles[2, , ] & truth <= ir$quantiles[4, , ]),
    coverage90 = mean(truth >= ir$quantiles[1, , ] & truth <= ir$quantiles[5, , ]),
    sign_probability_truth = mean(vapply(which(abs(truth) > 1e-10), function(j) {
      r <- (j - 1L) %% nrow(truth) + 1L; h <- (j - 1L) %/% nrow(truth) + 1L
      mean(ir$draws[r, h, ] * sign(truth[r, h]) > 0)
    }, numeric(1))),
    stable_rate = jcer_diagnostics(fit)$stable_rate, synthetic = TRUE)
}
acc <- do.call(rbind, results)
write.csv(acc, "benchmarks/results/accuracy_results.csv", row.names = FALSE)
for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf(paste0("benchmarks/results/accuracy.", ext), width = 8, height = 5) else png(paste0("benchmarks/results/accuracy.", ext), 1000, 650)
  barplot(colMeans(acc[, c("coverage68", "coverage90")]), ylim = c(0, 1), col = "steelblue",
    main = paste("SYNTHETIC DATA | conditional country inference |", nr, "replications"), ylab = "Empirical pointwise IRF coverage")
  abline(h = c(.68, .90), lty = 3); dev.off()
}
print(comparison); print(jcer_compare_irf(refir, newir)); print(colMeans(acc[, 2:7]))
