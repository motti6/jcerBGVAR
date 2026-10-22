Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
sim <- jcer_simulate_gvar(n_countries = 3, n_obs = 120, seed = 42, contemporaneous = FALSE)
fit <- jcer_bgvar(sim$data, sim$W, method = "conjugate", prior = "mniw",
                  lags = 1, draws = 1000, seed = 42, contemporaneous = FALSE)
print(summary(fit))
ir <- jcer_irf(fit, impulse = "US.r_short", response = c("JP.r_short", "JP.y"),
               identification = "generalized", horizon = 12)
dir.create("artifacts", showWarnings = FALSE)
pdf("artifacts/synthetic_spillovers.pdf", width = 8, height = 7); plot(ir); dev.off()
png("artifacts/synthetic_spillovers.png", width = 1000, height = 850); plot(ir); dev.off()
stopifnot(all(is.finite(ir$draws)))
cat("README example completed in a fresh R session.\n")
