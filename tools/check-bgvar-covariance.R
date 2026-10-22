Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(BGVAR)
b <- readRDS("artifacts/historical-pilot.rds")
S <- b$stacked.results$S_large[, , 1]; print(S[1:4, 1:4]); cat("Symmetry:", max(abs(S - t(S))), "\n")
for (fn in c(".irf.chol", ".irf.sign.zero")) {
  txt <- deparse(body(get(fn, asNamespace("BGVAR"))))
  hit <- grep("Smat|chol|Sigma|sigma|Ginv|impact|Pmat|S_c", txt)
  ix <- sort(unique(unlist(lapply(hit, function(j) seq(max(1,j-2), min(length(txt),j+3))))))
  cat("\n", fn, "\n", paste(txt[ix], collapse = "\n"), "\n")
}
