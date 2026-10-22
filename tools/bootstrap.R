Sys.setlocale("LC_CTYPE", ".UTF-8")
dir.create(".library", showWarnings = FALSE)
.libPaths(c(normalizePath(".library"), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 120)
for (p in c("roxygen2", "BGVAR")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    tryCatch(install.packages(p, lib = ".library", type = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType")), error = function(e) message(conditionMessage(e)))
  }
}
print(installed.packages()[intersect(c("BGVAR", "roxygen2"), rownames(installed.packages())), c("Package", "Version"), drop = FALSE])
