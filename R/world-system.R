.world_system <- function(B, Sigma, prepared, condition_tol = 1e-10) {
  m <- length(prepared$labels); k <- length(prepared$variables); p <- prepared$lags
  G0 <- matrix(0, m, m); G <- replicate(p, matrix(0, m, m), simplify = FALSE)
  Se <- matrix(0, m, m); cc <- numeric(m)
  for (i in seq_along(B)) {
    b <- B[[i]]; d <- prepared$blocks[[i]]; idx <- d$index
    cc[idx] <- b[1L, ]; G0[idx, ] <- d$E
    off <- 1L + k * p
    if (prepared$contemporaneous) {
      G0[idx, ] <- G0[idx, , drop = FALSE] - t(b[off + seq_len(k), , drop = FALSE]) %*% d$H
      off <- off + k
    }
    for (l in seq_len(p)) G[[l]][idx, ] <-
      t(b[1L + (l - 1L) * k + seq_len(k), , drop = FALSE]) %*% d$E +
      t(b[off + (l - 1L) * k + seq_len(k), , drop = FALSE]) %*% d$H
    .spd(Sigma[[i]], "Sigma_i"); Se[idx, idx] <- Sigma[[i]]
  }
  rc <- rcond(G0)
  if (!is.finite(rc) || rc < condition_tol) stop("G0 is singular or ill-conditioned")
  fac <- qr(G0, LAPACK = TRUE)
  A <- lapply(G, function(g) qr.coef(fac, g))
  # Transform a covariance factor, never explicitly invert G0.
  L <- qr.coef(fac, t(chol(Se)))
  Su <- tcrossprod(L)
  comp <- matrix(0, m * p, m * p); comp[seq_len(m), ] <- do.call(cbind, A)
  if (p > 1L) comp[m + seq_len(m * (p - 1L)), seq_len(m * (p - 1L))] <- diag(m * (p - 1L))
  rho <- max(Mod(eigen(comp, only.values = TRUE)$values))
  list(G0 = G0, G = G, Sigma_e = Se, Sigma_u = Su, c = cc, A = A,
       factor = fac, rcond = rc, spectral_radius = rho, stable = rho < 1)
}
.recur_irf <- function(world, impact, horizon) {
  m <- nrow(world$G0); out <- matrix(0, m, horizon + 1L); out[, 1L] <- impact
  if (horizon > 0L) for (h in seq_len(horizon)) {
    rhs <- numeric(m)
    for (l in seq_len(min(h, length(world$G)))) rhs <- rhs + world$G[[l]] %*% out[, h - l + 1L]
    out[, h + 1L] <- qr.coef(world$factor, rhs)
  }
  out
}
