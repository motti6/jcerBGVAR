Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
dir.create("benchmarks/results", recursive = TRUE, showWarnings = FALSE)
full <- "--full" %in% commandArgs(TRUE)
b <- jcer_benchmark(n_countries = c(3, 6, 12, 20), n_obs = c(80, 120, 200),
  lags = c(1, 2), draws = if (full) 200 else 100, repeats = 5,
  cores = if (full) c(1, 2, 4, 8) else 1, output_dir = "benchmarks/results")
if (!full) {
  scaling <- jcer_benchmark(n_countries = 12, n_obs = 120, lags = 1, draws = 100, repeats = 5, cores = c(1, 2, 4, 8))
  write.csv(scaling, "benchmarks/results/core_scaling.csv", row.names = FALSE)
}
writeLines(c(capture.output(sessionInfo()), paste("Logical cores:", parallel::detectCores()),
  paste("CPU:", Sys.getenv("PROCESSOR_IDENTIFIER")), paste("BLAS:", extSoftVersion()["BLAS"]),
  paste("OMP_NUM_THREADS:", Sys.getenv("OMP_NUM_THREADS", unset = "unset")),
  paste("OPENBLAS_NUM_THREADS:", Sys.getenv("OPENBLAS_NUM_THREADS", unset = "unset")),
  "Peak RAM not measured; fit_bytes is serialized R object footprint, not peak process RAM.",
  "Load/compile time excluded from warmed stage timings. Pure R: no package compilation."), "benchmarks/results/environment.txt")
print(aggregate(b[, c("country_estimation", "world_construction", "irf", "total")], b["n_countries"], median))
