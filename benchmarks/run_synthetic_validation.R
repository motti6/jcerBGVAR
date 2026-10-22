Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
source("benchmarks/lag-only-validation.R")
outdir <- "benchmarks/results/synthetic_revised_20261022"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pilot <- "--pilot" %in% commandArgs(TRUE)
scenarios <- data.frame(case = c("lag_only_moderate", "lag_only_persistent", "simultaneous_diagnostic"),
  persistence = c(.4, .75, .4), contemporaneous = c(FALSE, FALSE, TRUE),
  replications = if (pilot) rep(5L, 3) else c(500L, 500L, 200L),
  role = c("matched conditional/full lag-only likelihood", "matched conditional/full lag-only likelihood",
           "deliberate weak-exogeneity violation diagnostic; not a validity certificate"))
draws <- if (pilot) 50L else 250L
rows <- list(); coefficient <- list()
for (ci in seq_len(nrow(scenarios))) {
  cfg <- scenarios[ci, ]
  for (rep in seq_len(cfg$replications)) {
    seed <- 20000L + ci * 1000L + rep
    s <- jcer_simulate_gvar(3, 120, seed = seed, contemporaneous = cfg$contemporaneous,
                            persistence = cfg$persistence)
    if (cfg$contemporaneous) {
      f <- jcer_bgvar(s$data, s$W, draws = draws, seed = seed + 10000L,
                       contemporaneous = TRUE, stability = "keep")
      ir <- jcer_irf(f, "US.r_short", c("JP.r_short", "JP.y"), horizon = 8)
      stable_rate <- jcer_diagnostics(f)$stable_rate
    } else {
      ir <- lag_only_validation(s, draws, seed + 10000L)
      f <- list(countries = ir$countries); stable_rate <- ir$stable_rate
    }
    truth <- s$true_irf[match(ir$response, rownames(s$true_irf)), 1:9, drop = FALSE]
    # Report each response/horizon separately. h=0 is omitted from coverage:
    # lag-only foreign responses are mechanically zero on impact.
    for (r in seq_along(ir$response)) for (h in 1:8) {
      q <- ir$quantiles[, r, h + 1L]; tr <- truth[r, h + 1L]
      rows[[length(rows) + 1L]] <- data.frame(case = cfg$case, role = cfg$role, replication = rep,
        response = ir$response[r], horizon = h, truth = tr, median = q[3], error = q[3] - tr,
        cover68 = tr >= q[2] && tr <= q[4], cover90 = tr >= q[1] && tr <= q[5],
        width68 = q[4] - q[2], width90 = q[5] - q[1],
        sign_probability = if (abs(tr) > 1e-10) mean(ir$draws[r, h + 1L, ] * sign(tr) > 0) else NA_real_,
        stable_rate = stable_rate, posterior_draws = draws,
        contemporaneous = cfg$contemporaneous, persistence = cfg$persistence)
    }
    # Domestic own-lag rate coefficient coverage, independently of world assembly.
    br <- 1L + 2L; trb <- s$truth$B[[1]][br, 2L]
    q <- quantile(f$countries[[1]]$B[br, 2, ], c(.05, .16, .5, .84, .95), names = FALSE)
    coefficient[[length(coefficient) + 1L]] <- data.frame(case = cfg$case, replication = rep, truth = trb,
      median = q[3], cover68 = trb >= q[2] && trb <= q[4], cover90 = trb >= q[1] && trb <= q[5])
    if (rep %% 25L == 0L || rep == cfg$replications) {
      write.csv(do.call(rbind, rows), file.path(outdir, "pointwise_raw.csv"), row.names = FALSE)
      write.csv(do.call(rbind, coefficient), file.path(outdir, "coefficient_raw.csv"), row.names = FALSE)
      message(sprintf("%s %d/%d finished (%d independent datasets; %d posterior draws each)",
        cfg$case, rep, cfg$replications, rep, draws))
    }
  }
}
raw <- do.call(rbind, rows)
groups <- split(raw, interaction(raw$case, raw$response, raw$horizon, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(z) {
  n <- nrow(z); p68 <- mean(z$cover68); p90 <- mean(z$cover90)
  data.frame(case = z$case[1], role = z$role[1], response = z$response[1], horizon = z$horizon[1],
    independent_replications = n, bias = mean(z$error), bias_mcse = sd(z$error) / sqrt(n),
    rmse = sqrt(mean(z$error^2)), coverage68 = p68, coverage68_mcse = sqrt(p68 * (1 - p68) / n),
    coverage90 = p90, coverage90_mcse = sqrt(p90 * (1 - p90) / n),
    mean_width90 = mean(z$width90), mean_correct_sign_probability = mean(z$sign_probability, na.rm = TRUE),
    stable_rate = mean(z$stable_rate))
}))
write.csv(summary, file.path(outdir, "pointwise_summary.csv"), row.names = FALSE)
# Aggregate standard errors use independent dataset-level averages: horizons and
# responses in a dataset are correlated and are not counted as extra replications.
per_dataset <- aggregate(raw[, c("error", "cover68", "cover90", "stable_rate")], raw[, c("case", "replication")], mean)
agg <- do.call(rbind, lapply(split(per_dataset, per_dataset$case), function(z) data.frame(case = z$case[1],
  independent_replications = nrow(z), coverage68 = mean(z$cover68), coverage68_mcse = sd(z$cover68) / sqrt(nrow(z)),
  coverage90 = mean(z$cover90), coverage90_mcse = sd(z$cover90) / sqrt(nrow(z)),
  bias = mean(z$error), bias_mcse = sd(z$error) / sqrt(nrow(z)), stable_rate = mean(z$stable_rate))))
write.csv(agg, file.path(outdir, "scenario_summary.csv"), row.names = FALSE)
coef <- do.call(rbind, coefficient)
write.csv(aggregate(coef[, c("cover68", "cover90")], coef["case"], mean), file.path(outdir, "coefficient_summary.csv"), row.names = FALSE)
writeLines(c("Lag-only scenarios set B_i0=0 and omit current foreign regressors in estimation.",
  "With fixed initial lags and independent country errors, the joint lag-only likelihood factorizes by country.",
  "The simultaneous case deliberately retains the previous misspecification risk as a diagnostic, not validation.",
  "All figures/data are synthetic. h0 foreign responses excluded from coverage.",
  "Pointwise MCSE uses independent datasets; aggregate MCSE clusters response/horizon cells by dataset.",
  paste("Posterior draws per dataset:", draws), capture.output(sessionInfo())), file.path(outdir, "design.txt"))
pdf(file.path(outdir, "pointwise_coverage.pdf"), 9, 6)
par(mfrow = c(3, 1), mar = c(3, 4, 3, 1))
for (cs in scenarios$case) {
  z <- subset(summary, case == cs & response == "JP.r_short")
  plot(z$horizon, z$coverage90, ylim = c(0, 1), type = "b", xlab = "Horizon", ylab = "90% coverage",
       main = paste("SYNTHETIC DATA |", cs))
  arrows(z$horizon, pmax(0, z$coverage90 - 1.96 * z$coverage90_mcse),
         z$horizon, pmin(1, z$coverage90 + 1.96 * z$coverage90_mcse), angle = 90, code = 3, length = .04)
  abline(h = .90, lty = 3)
}
dev.off()
print(agg)
