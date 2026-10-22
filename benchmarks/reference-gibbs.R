# Independent implementation for validation only, deliberately using dense inverse
# formulas rather than production posterior helpers. Same MNIW conditional model.
reference_gibbs <- function(Y, Z, prior, draws = 2000, burnin = 1000, seed = 1) {
  set.seed(seed)
  k <- ncol(Y); z <- ncol(Z); n <- nrow(Y)
  Q0 <- solve(prior$V0); Qn <- Q0 + crossprod(Z); V <- solve(Qn)
  M <- V %*% (Q0 %*% prior$B0 + crossprod(Z, Y))
  S <- prior$S0 + crossprod(Y) + t(prior$B0) %*% Q0 %*% prior$B0 - t(M) %*% Qn %*% M
  Sigma <- diag(k); B <- M; out <- array(0, c(z, k, draws)); ss <- array(0, c(k, k, draws))
  for (i in seq_len(draws + burnin)) {
    B <- M + t(chol(V)) %*% matrix(rnorm(z * k), z, k) %*% chol(Sigma)
    scale <- S + crossprod(B - M, Qn %*% (B - M))
    Sigma <- solve(matrix(stats::rWishart(1, prior$nu0 + n + z, solve(scale))[, , 1], k, k))
    if (i > burnin) { out[, , i - burnin] <- B; ss[, , i - burnin] <- Sigma }
  }
  list(B = out, Sigma = ss)
}
chain_ess <- function(x) {
  aa <- as.numeric(stats::acf(x, plot = FALSE, lag.max = min(1000, length(x) - 1))$acf)[-1]
  pairs <- aa[seq(1, length(aa) - 1, 2)] + aa[seq(2, length(aa), 2)]
  first <- which(pairs <= 0)[1]; if (!is.na(first)) pairs <- if (first > 1) pairs[seq_len(first - 1)] else numeric()
  min(length(x), length(x) / max(1, 1 + 2 * sum(pairs)))
}
classic_rhat <- function(chains) {
  n <- nrow(chains); W <- mean(apply(chains, 2, var)); B <- n * var(colMeans(chains))
  sqrt(((n - 1) / n * W + B / n) / W)
}
