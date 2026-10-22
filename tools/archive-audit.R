Sys.setlocale("LC_CTYPE", ".UTF-8")
version <- read.dcf("DESCRIPTION")[1, "Version"]
archive <- paste0("artifacts/preJCERBGVAR_", version, ".tar.gz")
files <- utils::untar(archive, list = TRUE)
stopifnot(!any(grepl("/(\\.library|artifacts|benchmarks|data-raw|docs)/", files)),
          !any(grepl("eerData|W.trade0012|bgvar_body_inspection", files)))
stopifnot(sum(grepl("/inst/doc/.*\\.html$", files)) == 4)
writeLines(files, "artifacts/source_archive_contents.txt")
cat("Archive audited:", length(files), "entries; four rendered HTML vignettes; no local libraries or BGVAR data/code logs.\n")
