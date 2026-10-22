Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(BGVAR)
e <- new.env(); data("eerData", package = "BGVAR", envir = e)
cc <- c("US", "JP", "EA", "KR"); vv <- c("y", "Dp", "stir", "ltir")
d <- lapply(e$eerData[cc], function(x) {
  out <- cbind(y = 100 * diff(x[, "y"]), Dp = 100 * x[-1, "Dp"],
               stir = 100 * x[-1, "stir"], ltir = 100 * x[-1, "ltir"])
  rownames(out) <- rownames(x)[-1]; out
})
w <- e$W.trade0012[cc, cc]; w <- w / rowSums(w)
set.seed(1)
m <- bgvar(d, w, plag = 1, draws = 20, burnin = 20, prior = "MN", SV = FALSE,
           eigen = FALSE, expert = list(save.country.store = TRUE), verbose = FALSE)
str(m$stacked.results, max.level = 1)
str(m$cc.results$store[[1]], max.level = 1)
str(m$cc.results[c("coeffs", "sig")], max.level = 2)
sh <- get_shockinfo("girf"); sh$shock <- "US.stir"; sh$scale <- .25
g <- irf(m, n.ahead = 12, shockinfo = sh, quantiles = c(.05, .16, .5, .84, .95),
         expert = list(save.store = TRUE), verbose = FALSE)
print(dim(g$IRF_store)); print(dimnames(g$IRF_store)); print(g$IRF_store["US.stir", 1, 1, 1:5])
print(dim(g$posterior)); print(dimnames(g$posterior))
for (fn in c(".irf.girf", ".gvar.stacking.wrapper")) { cat("\n", fn, "\n"); print(get(fn, asNamespace("BGVAR"))) }
bodyir <- deparse(body(get("irf.bgvar", asNamespace("BGVAR"))))
cat(paste(bodyir[520:542], collapse = "\n"))
saveRDS(m, "artifacts/historical-pilot.rds")
saveRDS(g, "artifacts/historical-pilot-irf.rds")
