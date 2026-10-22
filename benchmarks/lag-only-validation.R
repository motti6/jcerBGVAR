# Specialized validation harness for known G0=I and p=1. Uses production MNIW
# country draws, but independently assembles reduced matrices. It is checked
# against the general fitter/IRF for shared seeds before Monte Carlo execution.
# No timing claims use this validation harness.
lag_only_validation <- function(s, draws, seed) {
  stopifnot(identical(s$contemporaneous, FALSE))
  p <- jcer_prepare_data(s$data, s$W, contemporaneous = FALSE)
  n <- length(p$countries); k <- length(p$variables); m <- n * k
  tasks <- preJCERBGVAR:::.with_seed(seed, lapply(p$blocks, function(z) list(
    Y = z$Y, Z = z$Z, prior = preJCERBGVAR:::.prior_mniw(ncol(z$Z), k, list()),
    seed = sample.int(.Machine$integer.max, 1), draws = draws)))
  cc <- lapply(tasks, preJCERBGVAR:::.country_draws)
  rr <- match(c("JP.r_short", "JP.y"), p$labels); jj <- match("US.r_short", p$labels)
  arr <- array(0, c(2, 9, draws), dimnames = list(c("JP.r_short", "JP.y"), as.character(0:8), NULL))
  stable <- logical(draws)
  for (d in seq_len(draws)) {
    F <- matrix(0, m, m); Se <- matrix(0, m, m)
    for (i in seq_len(n)) {
      idx <- (i - 1L) * k + seq_len(k); B <- cc[[i]]$B[, , d]
      F[idx, idx] <- t(B[1L + seq_len(k), , drop = FALSE])
      for (j in seq_len(n)) {
        jdx <- (j - 1L) * k + seq_len(k)
        F[idx, jdx] <- F[idx, jdx, drop = FALSE] + s$W[i, j] * t(B[1L + k + seq_len(k), , drop = FALSE])
      }
      Se[idx, idx] <- cc[[i]]$Sigma[, , d]
    }
    stable[d] <- max(Mod(eigen(F, only.values = TRUE)$values)) < 1
    path <- matrix(0, m, 9); path[, 1] <- Se[, jj] / sqrt(Se[jj, jj])
    for (h in 1:8) path[, h + 1L] <- F %*% path[, h]
    arr[, , d] <- path[rr, , drop = FALSE]
  }
  list(draws = arr, quantiles = apply(arr, c(1, 2), quantile, probs = c(.05, .16, .5, .84, .95), names = FALSE),
       response = dimnames(arr)[[1]], stable_rate = mean(stable), countries = cc)
}
lag_only_validation <- compiler::cmpfun(lag_only_validation)
