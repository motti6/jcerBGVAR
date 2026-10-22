#' Benchmark conditional GVAR estimation and selected IRFs
#' @param n_countries Vector of country counts.
#' @param n_obs Vector of sample sizes.
#' @param lags Vector of lag orders.
#' @param draws Independent draws per fit.
#' @param repeats Timed repetitions, following one warm-up per configuration.
#' @param cores Vector of country worker counts.
#' @param seed Base seed.
#' @param output_dir Optional directory for raw CSV and plots.
#' @return Data frame with measured stage times, fit size and stability. Peak RAM is NA unless independently measured.
#' @export
#' @examples
#' b <- jcer_benchmark(n_countries = 3, n_obs = 40, draws = 5, repeats = 1)
jcer_benchmark <- function(n_countries = c(3, 6, 12, 20), n_obs = c(80, 120, 200),
                            lags = c(1, 2), draws = 100, repeats = 5, cores = 1,
                            seed = 42, output_dir = NULL) {
  repeats <- .integer_arg(repeats, "repeats")
  grid <- expand.grid(n_countries = n_countries, n_obs = n_obs, lags = lags, cores = cores)
  rows <- list(); count <- 0L
  for (g in seq_len(nrow(grid))) {
    cfg <- grid[g, ]; s <- jcer_simulate_gvar(cfg$n_countries, cfg$n_obs, seed, lags = cfg$lags)
    invisible(jcer_bgvar(s$data, s$W, lags = cfg$lags, draws = 2, seed = seed, cores = cfg$cores))
    for (r in seq_len(repeats)) {
      f <- jcer_bgvar(s$data, s$W, lags = cfg$lags, draws = draws, seed = seed + r, cores = cfg$cores)
      t0 <- proc.time()[[3L]]
      for (w in f$worlds) invisible(w$Sigma_u[, 2] / sqrt(w$Sigma_u[2, 2]))
      identification <- proc.time()[[3L]] - t0
      irf_time <- system.time(ir <- jcer_irf(f, "US.r_short", c("JP.r_short", "JP.y"), horizon = 12))[[3L]]
      tmp <- tempfile(fileext = ".pdf"); grDevices::pdf(tmp)
      plot_time <- system.time(plot(ir))[[3L]]; grDevices::dev.off(); unlink(tmp)
      count <- count + 1L
      rows[[count]] <- cbind(cfg, data.frame(repetition = r, draws = draws, preprocessing = unname(f$timing[1]),
        country_estimation = unname(f$timing[2]), world_construction = unname(f$timing[3]),
        identification = identification, irf = irf_time, visualization = plot_time,
        total = unname(f$timing[4]) + irf_time + plot_time, fit_bytes = as.numeric(object.size(f)), peak_ram_bytes = NA_real_,
        stable_rate = jcer_diagnostics(f)$stable_rate, synthetic = TRUE))
    }
  }
  ans <- do.call(rbind, rows)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(ans, file.path(output_dir, "benchmark_results.csv"), row.names = FALSE)
    for (ext in c("pdf", "png")) {
      path <- file.path(output_dir, paste0("benchmark_stages.", ext))
      if (ext == "pdf") grDevices::pdf(path, width = 9, height = 5) else grDevices::png(path, width = 1200, height = 700)
      med <- stats::aggregate(ans[, c("country_estimation", "world_construction", "irf")], ans[, c("n_countries", "n_obs", "lags", "cores")], stats::median)
      graphics::barplot(t(as.matrix(med[, c("country_estimation", "world_construction", "irf")])), beside = FALSE,
        names.arg = paste(med$n_countries, med$n_obs, med$lags, med$cores, sep = "/"), col = c("steelblue", "orange", "grey50"),
        ylab = "Median elapsed seconds", main = "SYNTHETIC DATA | country count / T / lags / cores")
      graphics::legend("topleft", c("Country draws", "World construction", "Selected IRF"), fill = c("steelblue", "orange", "grey50"))
      grDevices::dev.off()
    }
  }
  ans
}
