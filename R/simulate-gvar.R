#' Simulate a stable GVAR with known world dynamics
#' @param n_countries Number of countries, at least two (US, JP, EA, KR, then C5, ...).
#' @param n_obs Retained observations.
#' @param seed Random seed.
#' @param variables Common variable names. Defaults to GDP and short rate analogues.
#' @param lags True lag order.
#' @param rate_unit Rate units; decimal output rescales all r_short/r_long quantities.
#' @param contemporaneous Include simultaneous foreign transmission. Set FALSE
#'   for a lag-only DGP consistent with conditional regression exogeneity.
#' @param persistence Domestic diagonal persistence, divided equally across lags.
#' @param foreign_lag Foreign diagonal lag coefficient, divided across lags.
#' @return List with synthetic country data, W, truth (world system and coefficient blocks), true generalized IRF and units.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(seed = 42)
#' s$truth$world$spectral_radius
jcer_simulate_gvar <- function(n_countries = 3, n_obs = 120, seed = 1,
                               variables = c("y", "r_short"), lags = 1,
                               rate_unit = c("percent", "decimal"), contemporaneous = TRUE,
                               persistence = .4, foreign_lag = .05) {
  n <- .integer_arg(n_countries, "n_countries", 2L); n_obs <- .integer_arg(n_obs, "n_obs", 10L)
  lags <- .integer_arg(lags, "lags"); rate_unit <- match.arg(rate_unit)
  if (!is.logical(contemporaneous) || length(contemporaneous) != 1 || is.na(contemporaneous)) stop("contemporaneous must be TRUE or FALSE")
  if (length(persistence) != 1 || !is.finite(persistence) || length(foreign_lag) != 1 || !is.finite(foreign_lag)) stop("finite scalar dynamic coefficients required")
  if (!.named_unique(variables) || !length(variables) || any(grepl("\\.", variables))) stop("variables must be unique names without dots")
  countries <- c("US", "JP", "EA", "KR", if (n > 4) paste0("C", seq.int(5, n)))[seq_len(n)]
  k <- length(variables); m <- n * k
  W <- matrix(1 / (n - 1), n, n, dimnames = list(countries, countries)); diag(W) <- 0
  blank <- setNames(lapply(seq_len(n), function(i) matrix(0, n_obs, k, dimnames = list(NULL, variables))), countries)
  pre <- jcer_prepare_data(blank, W, lags, contemporaneous = contemporaneous, rate_unit = rate_unit)
  z <- ncol(pre$blocks[[1]]$Z); B <- rep(list(matrix(0, z, k)), n)
  scale <- rep(1, k); if (rate_unit == "decimal") scale[variables %in% c("r_short", "r_long")] <- 0.01
  D <- diag(scale, k); Di <- diag(1 / scale, k)
  for (i in seq_len(n)) {
    for (l in seq_len(lags)) {
      A <- diag(persistence / lags, k)
      if (all(c("y", "r_short") %in% variables)) A[match("y", variables), match("r_short", variables)] <- -0.06 / lags
      B[[i]][1L + (l - 1L) * k + seq_len(k), ] <- t(D %*% A %*% Di)
      B[[i]][1L + k * lags + (if (contemporaneous) k else 0L) + (l - 1L) * k + seq_len(k), ] <- diag(foreign_lag / lags, k)
    }
    if (contemporaneous) B[[i]][1L + k * lags + seq_len(k), ] <- diag(0.08, k)
  }
  S <- rep(list(diag((0.2 * scale)^2, k)), n)
  w <- .world_system(B, S, pre)
  if (!w$stable) stop("requested synthetic dynamics are not stable")
  X <- .with_seed(seed, {
    x <- matrix(0, n_obs + 200L + lags, m)
    for (t in seq.int(lags + 1L, nrow(x))) {
      rhs <- numeric(m)
      for (l in seq_len(lags)) rhs <- rhs + w$G[[l]] %*% x[t - l, ]
      x[t, ] <- qr.coef(w$factor, rhs + drop(t(chol(w$Sigma_e)) %*% stats::rnorm(m)))
    }
    tail(x, n_obs)
  })
  data <- setNames(lapply(seq_len(n), function(i) {
    out <- X[, (i - 1L) * k + seq_len(k), drop = FALSE]; colnames(out) <- variables; out
  }), countries)
  attr(data, "synthetic") <- TRUE
  j <- if ("r_short" %in% variables) match("r_short", variables) else 1L
  ir <- .recur_irf(w, w$Sigma_u[, j] / sqrt(w$Sigma_u[j, j]), 20)
  rownames(ir) <- pre$labels
  list(data = data, W = W, truth = list(B = B, Sigma = S, world = w),
       true_irf = ir, rate_unit = rate_unit, synthetic = TRUE,
       contemporaneous = contemporaneous,
       exogeneity = if (contemporaneous) "simultaneous DGP: conditional foreign exogeneity is not guaranteed" else
         "lag-only DGP: innovations independent of all lagged domestic and foreign regressors")
}
