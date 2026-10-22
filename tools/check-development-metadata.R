Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
d <- utils::packageDescription("preJCERBGVAR")
stopifnot(d$Package == "preJCERBGVAR", d$Version == "0.1.1.9000", d$License == "GPL-3",
          grepl("motohiro watanabe", d$Maintainer, fixed = TRUE),
          grepl("watanabemotosan@gmail.com", d$Maintainer, fixed = TRUE))
print(d[c("Package", "Version", "Author", "Maintainer", "License", "URL")])
print(utils::citation("preJCERBGVAR"))
cat("Development metadata and installed R citation validated.\n")
