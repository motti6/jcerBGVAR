Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
if (requireNamespace("RcppParallel", quietly = TRUE)) RcppParallel::setThreadOptions(numThreads = 1)
source("benchmarks/historical-common.R")
outdir <- "benchmarks/results/historical_20261022"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(TRUE)
pilot <- "--pilot" %in% args
draws <- if (pilot) 50L else 1000L; burnin <- if (pilot) 50L else 1000L
repeats <- if (pilot) 1L else 5L
configs <- data.frame(case = c("NG_4", "NG_12", "NG_20", "MN_4_long", "NG_4_SV"),
  n = c(4L, 12L, 20L, 4L, 4L), long = c(FALSE, FALSE, FALSE, TRUE, FALSE),
  prior = c("NG", "NG", "NG", "MN", "NG"), SV = c(FALSE, FALSE, FALSE, FALSE, TRUE))
if (pilot) configs <- configs[c(1, 4), ]
rows <- list(); model_rows <- list(); definition_rows <- list(); convergence <- list()
flush <- function() {
  if (length(rows)) write.csv(do.call(rbind, rows), file.path(outdir, "same_bgvar_timing_raw.csv"), row.names = FALSE)
  if (length(model_rows)) write.csv(do.call(rbind, model_rows), file.path(outdir, "different_model_timing_raw.csv"), row.names = FALSE)
  if (length(definition_rows)) write.csv(do.call(rbind, definition_rows), file.path(outdir, "irf_distribution_differences.csv"), row.names = FALSE)
  if (length(convergence)) write.csv(do.call(rbind, convergence), file.path(outdir, "selected_draw_diagnostics.csv"), row.names = FALSE)
}
source("benchmarks/reference-gibbs.R")
for (ci in seq_len(nrow(configs))) {
  cfg <- configs[ci, ]; d <- historical_data(cfg$n, cfg$long)
  manifest <- data.frame(country = d$countries, selected_trade_mass = d$retained_trade_mass,
    excluded_trade_mass = 1 - d$retained_trade_mass, period_start = d$period[1], period_end = d$period[2],
    variables = paste(d$variables, collapse = ";"), transform = "100*diff(logGDP);100*Dp;100*decimal_market_rates",
    explicit_input_rate_unit = "decimal", output_rate_unit = "percent", n_obs = nrow(d$data[[1]]))
  write.csv(manifest, file.path(outdir, paste0(cfg$case, "_manifest.csv")), row.names = FALSE)
  # Warm-up runs are excluded; loading and data transformations are excluded for both methods.
  invisible(capture.output(warm <- historical_bgvar_fit(d, min(50L, draws), min(50L, burnin), 9, cfg$prior, cfg$SV)))
  message("Warm-up ", cfg$case, ": ", paste(colnames(warm$xglobal), collapse = ","))
  invisible(capture.output(historical_bgvar_irf(warm)))
  invisible(jcer_irf_bgvar(warm, "US.stir", horizon = 12))
  for (r in seq_len(repeats)) {
    seed <- 7300L + 100L * ci + r
    tf <- system.time(invisible(capture.output(fit <- historical_bgvar_fit(d, draws, burnin, seed, cfg$prior, cfg$SV))))[[3L]]
    # Alternate order to reduce thermal/cache/order bias. Both return all responses,
    # h=0..12, one shock, raw draws and the same five posterior quantiles.
    if (r %% 2L) {
      tb <- system.time(invisible(capture.output(ref <- historical_bgvar_irf(fit))))[[3L]]
      tn <- system.time(fast <- jcer_irf_bgvar(fit, "US.stir", horizon = 12))[[3L]]
    } else {
      tn <- system.time(fast <- jcer_irf_bgvar(fit, "US.stir", horizon = 12))[[3L]]
      tb <- system.time(invisible(capture.output(ref <- historical_bgvar_irf(fit))))[[3L]]
    }
    raw <- historical_irf_array(ref)
    delta <- max(abs(raw - fast$draws)); relative <- delta / max(1, max(abs(raw)))
    qref <- aperm(array(ref$posterior[, , 1, ], c(dim(raw)[1:2], 5)), c(3, 1, 2))
    qdelta <- max(abs(qref - fast$quantiles))
    if (relative > 1e-10 || qdelta / max(1, max(abs(qref))) > 1e-10) stop("same-BGVAR IRF mismatch")
    rows[[length(rows) + 1L]] <- data.frame(case = cfg$case, repetition = r, countries = cfg$n,
      variables_per_country = length(d$variables), n_obs = nrow(d$data[[1]]), lags = 1, prior = cfg$prior, SV = cfg$SV,
      draws = draws, burnin = burnin, seed = seed, fit_seconds = tf, bgvar_irf_seconds = tb, fast_irf_seconds = tn,
      bgvar_total_seconds = tf + tb, fast_total_seconds = tf + tn,
      irf_speedup = tb / tn, end_to_end_speedup = (tf + tb) / (tf + tn),
      max_draw_absolute_error = delta, max_draw_relative_error = relative, max_quantile_absolute_error = qdelta,
      stable_draw_rate = mean(fit$stacked.results$F.eigen < 1),
      bgvar_result_bytes = as.numeric(object.size(ref)), fast_result_bytes = as.numeric(object.size(fast)),
      peak_ram_bytes = NA_real_, definition = "BGVAR country-error GIRF, unchanged posterior; 25bp own impact each draw")
    if (r == 1L) {
      # Separate shock-definition difference using identical posterior draws.
      world <- jcer_irf_bgvar(fit, "US.stir", horizon = 12, definition = "world_generalized")
      def_delta <- fast$quantiles - world$quantiles
      definition_rows[[length(definition_rows) + 1L]] <- data.frame(case = cfg$case,
        comparison = "same posterior, country-error GIRF vs world-generalized innovation: NOT numerical error",
        median_rmse = sqrt(mean(def_delta[3, , ]^2)), endpoint_mae = mean(abs(def_delta[c(1, 5), , ])))
      stable <- tryCatch(stable_bgvar_subset(fit), error = identity)
      if (!inherits(stable, "error")) {
        sref <- historical_bgvar_irf(stable); sfast <- jcer_irf_bgvar(stable, "US.stir", horizon = 12)
        serr <- max(abs(historical_irf_array(sref) - sfast$draws))
        write.csv(data.frame(case = cfg$case, requested_draws = draws, stable_draws = length(stable$stacked.results$F.eigen),
          max_absolute_error = serr, target = "both methods conditioned on strictly stable BGVAR world draws"),
          file.path(outdir, paste0(cfg$case, "_stable_sensitivity.csv")), row.names = FALSE)
      }
    }
    # Fit with different MNIW prior as a distinct practical model comparison.
    tm <- system.time(mniw <- jcer_bgvar(d$data, d$W, draws = draws, seed = seed,
      rate_unit = "percent", stability = "keep"))[[3L]]
    # Align world-generalized shock definitions before comparing model distributions.
    ti <- system.time(mir <- jcer_irf(mniw, "US.stir", fast$response, horizon = 12, shock_bp = 25))[[3L]]
    world <- jcer_irf_bgvar(fit, "US.stir", horizon = 12, definition = "world_generalized")
    diff <- mir$quantiles - world$quantiles
    model_rows[[length(model_rows) + 1L]] <- data.frame(case = cfg$case, repetition = r, countries = cfg$n,
      reference_prior = cfg$prior, reference_SV = cfg$SV, candidate_prior = "MNIW", candidate_SV = FALSE,
      draws = draws, bgvar_fit_seconds = tf, bgvar_world_irf_seconds = world$elapsed,
      mniw_fit_seconds = tm, mniw_irf_seconds = ti, mniw_total_seconds = tm + ti,
      practical_time_ratio = (tf + world$elapsed) / (tm + ti),
      median_irf_rmse = sqrt(mean(diff[3, , ]^2)), interval_endpoint_mae = mean(abs(diff[c(1, 5), , ])),
      median_sign_agreement = mean(sign(mir$quantiles[3, , ]) == sign(world$quantiles[3, , ])),
      mniw_stable_rate = jcer_diagnostics(mniw)$stable_rate,
      note = "prior/covariance model differs; discrepancy is not numerical approximation error")
    # Selected functional ESS is descriptive; speed at equal draw count does not imply equal MCMC precision.
    for (response in intersect(c("JP.y", "JP.Dp", "JP.stir", "JP.ltir"), fast$response)) {
      vals <- fast$draws[response, "4", ]
      convergence[[length(convergence) + 1L]] <- data.frame(case = cfg$case, repetition = r, functional = paste0(response, ":h4"),
        draws = draws, ESS = chain_ess(vals), note = "single-chain functional ESS only; not a convergence certificate")
    }
    if (r == repeats) {
      saveRDS(list(bgvar_same_definition = fast, bgvar_world_definition = world, mniw_world_definition = mir),
              file.path(outdir, paste0(cfg$case, "_irf_results.rds")))
      # No historical observations or complete BGVAR fits are exported.
    }
    flush()
    message(sprintf("%s repetition %d/%d: fit %.3fs; BGVAR IRF %.3fs; vector %.3fs; total ratio %.3f; error %.3g",
                    cfg$case, r, repeats, tf, tb, tn, (tf + tb) / (tf + tn), delta))
    rm(fit, ref, fast, mniw, mir, world); invisible(gc())
  }
}
writeLines(c(capture.output(sessionInfo()), paste("BGVAR:", packageVersion("BGVAR")),
  paste("Threads:", Sys.getenv("OMP_NUM_THREADS"), Sys.getenv("OPENBLAS_NUM_THREADS")),
  "All comparisons sequential; input data preparation and package loading excluded; warm-up excluded.",
  "Same-BGVAR timing includes extraction, shock construction, normalization, propagation and quantiles.",
  "BGVAR public IRF also produces auxiliary median structural summaries not returned by selected-vector IRF.",
  "No equal-model inference speedup claimed: inference is unchanged in the primary comparison.",
  "Historical data loaded from installed BGVAR; no observations copied into this repository."), file.path(outdir, "environment.txt"))
