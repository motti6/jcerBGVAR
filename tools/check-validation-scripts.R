Sys.setlocale("LC_CTYPE", ".UTF-8")
for (directory in c("R", "benchmarks", "tools")) {
  files <- list.files(directory, pattern = "\\.R$", full.names = TRUE)
  for (f in files) parse(file = f, encoding = "UTF-8")
  cat("Parsed", length(files), "R scripts in", directory, "\n")
}
