Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(BGVAR)
dir.create("artifacts/logs", recursive = TRUE, showWarnings = FALSE)
writeLines(c(capture.output(packageVersion("BGVAR")), capture.output(args(BGVAR::bgvar)),
  capture.output(args(BGVAR::irf)), capture.output(args(BGVAR::get_shockinfo)),
  capture.output(args(BGVAR::add_shockinfo))), "artifacts/logs/bgvar_api.txt")
writeLines(deparse(body(BGVAR::bgvar)), "artifacts/logs/bgvar_body_inspection.txt")
writeLines(deparse(body(getS3method("irf", "bgvar"))), "artifacts/logs/bgvar_irf_inspection.txt")
cat(readLines("artifacts/logs/bgvar_api.txt"), sep = "\n")
