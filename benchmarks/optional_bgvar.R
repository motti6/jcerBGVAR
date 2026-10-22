# Public BGVAR API adapter. Only this optional development benchmark uses tracing
# to measure a verified internal stacking boundary; no foreign source is copied.
Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
if (!requireNamespace("BGVAR", quietly = TRUE)) stop("Optional BGVAR package is not installed")
dir.create("benchmarks/results", recursive = TRUE, showWarnings = FALSE)
stopifnot(all(c("Data", "W", "plag", "draws", "burnin", "SV", "expert") %in% names(formals(BGVAR::bgvar))))
# Historical data loaded by reference, never stored in this project.
historical <- new.env(); utils::data("eerData", package = "BGVAR", envir = historical)
writeLines(c("Loaded BGVAR eerData and W.trade0012 in memory only.",
  paste("Available selected countries:", paste(intersect(c("US", "JP", "EA", "KR"), names(historical$eerData)), collapse = ", ")),
  "Historical 1995Q1-2013Q4 sample; stir is a short-term market rate, not the policy decision.",
  "No real-data causal IRF is claimed. Main benchmarks use identical synthetic observations."),
  "benchmarks/results/historical_data_availability.txt")
bench_bgvar <- function(draws = 100, repeats = 5) {
  s <- jcer_simulate_gvar(seed = 42)
  ns <- asNamespace("BGVAR")
  if (!exists(".gvar.stacking.wrapper", envir = ns, inherits = FALSE)) stop("BGVAR stacking boundary changed: re-inspect source")
  # Timers live only in this standalone benchmark process.
  .jcer_stage <- new.env(parent = emptyenv()); .jcer_stage$stack <- NA_real_
  assign(".jcer_stage", .jcer_stage, envir = .GlobalEnv)
  trace(".gvar.stacking.wrapper", where = ns,
    tracer = quote(.GlobalEnv$.jcer_stage$t0 <- proc.time()[[3L]]),
    exit = quote(.GlobalEnv$.jcer_stage$stack <- proc.time()[[3L]] - .GlobalEnv$.jcer_stage$t0), print = FALSE)
  on.exit({ untrace(".gvar.stacking.wrapper", where = ns); rm(".jcer_stage", envir = .GlobalEnv) })
  timed_apply <- function(X, FUN, ...) {
    t0 <- proc.time()[[3L]]; ans <- lapply(X, FUN, ...)
    .jcer_stage$country <- proc.time()[[3L]] - t0; ans
  }
  rows <- list()
  for (r in 0:repeats) {
    set.seed(42 + r)
    tt <- system.time(m <- BGVAR::bgvar(Data = s$data, W = s$W, plag = 1,
      draws = draws, burnin = 100, prior = "MN", SV = FALSE, thin = 1,
      eigen = FALSE, trend = FALSE, expert = list(applyfun = timed_apply), verbose = FALSE))[[3L]]
    shock <- BGVAR::get_shockinfo("girf"); shock$shock <- "US.r_short"; shock$scale <- .25
    gt <- system.time(g <- BGVAR::irf(m, n.ahead = 12, shockinfo = shock, verbose = FALSE))[[3L]]
    signs <- BGVAR::get_shockinfo("sign")
    signs <- BGVAR::add_shockinfo(signs, shock = "US.r_short", restriction = "US.r_short", sign = ">", horizon = 1, scale = .25, prob = 1)
    st <- system.time(si <- BGVAR::irf(m, n.ahead = 12, shockinfo = signs,
      expert = list(MaxTries = 1000), verbose = FALSE))[[3L]]
    if (r > 0) rows[[r]] <- data.frame(repetition = r, model = "BGVAR MN; SV=FALSE; model/prior differs from MNIW",
      version = as.character(utils::packageVersion("BGVAR")), draws = draws, burnin = 100,
      country_estimation = .jcer_stage$country, world_stacking = .jcer_stage$stack,
      other_fit = max(0, tt - .jcer_stage$country - .jcer_stage$stack), fit_total = tt,
      generalized_irf = gt, sign_identification_plus_irf = st,
      sign_only = NA_real_, synthetic = TRUE,
      notes = "Sign identification and IRF fused by public API; sign-only time unavailable, not estimated")
    if (r == repeats) saveRDS(list(generalized_summary = g$posterior, shockinfo = g$shockinfo,
      sign_rotation_count = si$rot.nr), "benchmarks/results/bgvar_synthetic_summary.rds")
  }
  ans <- do.call(rbind, rows)
  write.csv(ans, "benchmarks/results/bgvar_comparison.csv", row.names = FALSE)
  print(ans)
  ans
}
bench_bgvar()
