Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
source("benchmarks/lag-only-validation.R")
for (seed in c(3, 17, 41)) {
  s <- jcer_simulate_gvar(3, 120, seed, contemporaneous = FALSE, persistence = .75)
  f <- jcer_bgvar(s$data, s$W, draws = 100, seed = seed, contemporaneous = FALSE)
  a <- jcer_irf(f, "US.r_short", c("JP.r_short", "JP.y"), horizon = 8)
  b <- lag_only_validation(s, 100, seed)
  stopifnot(max(abs(a$draws - b$draws)) < 1e-12,
    max(abs(a$quantiles - b$quantiles)) < 1e-12,
    identical(f$countries, setNames(b$countries, names(f$countries))),
    jcer_diagnostics(f)$stable_rate == b$stable_rate)
}
cat("Independent G0=I validation assembly matches general fitter for three seeds and all 300 world draws.\n")
