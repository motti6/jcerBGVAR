#' Selected-vector IRFs from an existing BGVAR posterior
#'
#' Reuses an already fitted BGVAR object without changing its prior, country
#' sampler, covariance draws or stability selection. The default definition
#' matches BGVAR 2.7.0's country-error GIRF, whose impact is proportional to
#' Ginv times Sigma_e times e_j. This differs from a generalized innovation of
#' the reduced world system, which uses Sigma_u times e_j. Both are normalized
#' draw by draw to the specified own-variable impact.
#' @param model An installed BGVAR package's fitted object of class bgvar.
#' @param impulse One country.variable label.
#' @param response Response labels; NULL returns every variable in model order.
#' @param horizon Final horizon, including impact at zero.
#' @param shock_size Positive own-variable impact in the input unit. Default .25.
#' @param rate_unit Explicit rate unit, percent or decimal; no rescaling is done.
#' @param definition bgvar matches the fitted package's country-error GIRF;
#'   world_generalized uses the reduced-world innovation covariance.
#' @param normalization_tol Minimum absolute own impact before normalization.
#' @return JCERBGVAR_irf with draw arrays, five quantiles, definition, BGVAR
#'   version/retained-draw metadata and elapsed time. Setup, normalization,
#'   recurrence and quantile computation are included in elapsed time.
#' @export
#' @examples
#' if (requireNamespace("BGVAR", quietly = TRUE)) {
#'   s <- jcer_simulate_gvar(n_obs = 40)
#'   b <- BGVAR::bgvar(s$data, s$W, draws = 10, burnin = 10,
#'                     prior = "MN", SV = FALSE, eigen = FALSE, verbose = FALSE)
#'   ir <- jcer_irf_bgvar(b, "US.r_short", "JP.y", horizon = 3)
#' }
jcer_irf_bgvar <- function(model, impulse, response = NULL, horizon = 20,
                            shock_size = .25, rate_unit = c("percent", "decimal"),
                            definition = c("bgvar", "world_generalized"), normalization_tol = 1e-8) {
  start <- proc.time()[[3L]]
  if (!inherits(model, "bgvar")) stop("model must be a fitted BGVAR::bgvar object")
  rate_unit <- match.arg(rate_unit); definition <- match.arg(definition)
  st <- model$stacked.results; labels <- colnames(model$xglobal)
  m <- length(labels); p <- max(model$args$lags); nd <- length(st$F.eigen)
  if (is.null(response)) response <- labels
  proxy <- structure(list(prepared = list(labels = labels, rate_unit = rate_unit, synthetic = FALSE),
    method = "unchanged BGVAR posterior / selected-vector IRF",
    posterior_target = "BGVAR fitted country posterior and retained world draws, unchanged by post-processing"), class = "JCERBGVAR")
  a <- .irf_args(proxy, impulse, response, horizon, shock_size, normalization_tol)
  expected <- list(F_large = c(m, m, p, nd), S_large = c(m, m, nd), Ginv_large = c(m, m, nd))
  for (key in names(expected)) if (!identical(dim(st[[key]]), as.integer(expected[[key]])) || any(!is.finite(st[[key]])))
    stop("unsupported or nonfinite BGVAR stacked.results layout: ", key)
  values <- array(0, c(length(response), a$horizon + 1L, nd),
                  dimnames = list(response, as.character(0:a$horizon), NULL))
  for (d in seq_len(nd)) {
    gi <- st$Ginv_large[, , d]; se <- st$S_large[, , d]
    if (definition == "bgvar") impact <- drop(gi %*% se[, a$j]) else
      impact <- drop(gi %*% (se %*% gi[a$j, ]))
    if (abs(impact[a$j]) < normalization_tol) stop("near-zero BGVAR own impact; cannot normalize")
    impact <- impact * shock_size / impact[a$j]
    F <- lapply(seq_len(p), function(l) st$F_large[, , l, d])
    path <- .recur_reduced(F, impact, a$horizon)
    values[, , d] <- path[a$r, , drop = FALSE]
  }
  out <- .irf_result(values, proxy, impulse, response,
    if (definition == "bgvar") "BGVAR country-error GIRF" else "world generalized innovation",
    shock_size, list(definition = definition, retained_draws = nd,
      source_BGVAR_version = if (requireNamespace("BGVAR", quietly = TRUE)) as.character(utils::packageVersion("BGVAR")) else NA_character_,
      stable_draw_rate = mean(st$F.eigen < 1)))
  out$elapsed <- proc.time()[[3L]] - start
  out
}

.recur_reduced <- function(F, impact, horizon) {
  out <- matrix(0, length(impact), horizon + 1L); out[, 1L] <- impact
  if (horizon > 0L) for (h in seq_len(horizon)) for (l in seq_len(min(h, length(F))))
    out[, h + 1L] <- out[, h + 1L] + F[[l]] %*% out[, h - l + 1L]
  out
}
