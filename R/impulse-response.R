.irf_args <- function(model, impulse, response, horizon, shock_size, normalization_tol) {
  if (!inherits(model, "JCERBGVAR")) stop("model must be a JCERBGVAR fit")
  labels <- model$prepared$labels
  if (length(impulse) != 1L || !impulse %in% labels || !length(response) || any(!response %in% labels) || anyDuplicated(response)) stop("impulse/response must be unique country.variable labels in the model")
  horizon <- .integer_arg(horizon, "horizon", 0L)
  if (!is.null(shock_size) && (length(shock_size) != 1 || !is.numeric(shock_size) || !is.finite(shock_size) || shock_size <= 0)) stop("shock_size must be a positive finite impact in the input unit")
  if (length(normalization_tol) != 1 || !is.finite(normalization_tol) || normalization_tol <= 0) stop("normalization_tol must be positive")
  list(j = match(impulse, labels), r = match(response, labels), horizon = horizon)
}
.irf_result <- function(values, model, impulse, response, identification, shock_size, extra = list()) {
  quant <- apply(values, c(1, 2), stats::quantile, probs = c(.05, .16, .5, .84, .95), names = FALSE)
  structure(c(list(draws = values, quantiles = quant, impulse = impulse, response = response,
     horizon = seq_len(dim(values)[2]) - 1L, identification = identification,
     rate_unit = model$prepared$rate_unit, shock_size = shock_size, method = model$method,
     posterior_target = model$posterior_target, synthetic = model$prepared$synthetic,
     interpretation = if (identification == "generalized") "generalized innovation response; not an identified Fed policy causal effect" else
       "response conditional on the stated structural identification assumptions"), extra), class = "JCERBGVAR_irf")
}
#' Compute selected generalized or recursive impulse responses
#' @param model A JCERBGVAR fit.
#' @param impulse One country.variable innovation label.
#' @param response Vector of country.variable response labels.
#' @param horizon Last horizon, with impact at zero.
#' @param identification Either generalized or recursive.
#' @param shock_size Desired impact in the input unit, or NULL for one standard deviation.
#' @param shock_bp Optional positive basis point size, for an explicitly selected rate impulse.
#' @param ordering Full variable ordering for recursive identification; required for recursive IRFs.
#' @param normalization_tol Reject impacts below this absolute threshold before normalization.
#' @return JCERBGVAR_irf with response x horizon x draw array, quantiles and identification/unit metadata.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 40)
#' f <- jcer_bgvar(s$data, s$W, draws = 20, seed = 1)
#' ir <- jcer_irf(f, "US.r_short", c("JP.r_short", "JP.y"), horizon = 5, shock_bp = 25)
#' plot(ir)
jcer_irf <- function(model, impulse, response, horizon = 20,
                      identification = c("generalized", "recursive"), shock_size = NULL,
                      shock_bp = NULL, ordering = NULL, normalization_tol = 1e-8) {
  identification <- match.arg(identification)
  if (!is.null(shock_bp)) {
    if (!is.null(shock_size) || length(shock_bp) != 1 || !is.finite(shock_bp) || shock_bp <= 0) stop("use a positive shock_bp or shock_size, not both")
    if (!grepl("\\.(r_short|r_long|stir|ltir)$", impulse)) stop("shock_bp requires a rate impulse label")
    shock_size <- shock_bp / if (model$prepared$rate_unit == "percent") 100 else 10000
  }
  a <- .irf_args(model, impulse, response, horizon, shock_size, normalization_tol)
  labels <- model$prepared$labels; m <- length(labels)
  ord <- seq_len(m)
  if (identification == "recursive") {
    if (is.null(ordering) || length(ordering) != m || anyDuplicated(ordering) || !setequal(ordering, labels)) stop("recursive IRF requires an explicit complete ordering")
    ord <- match(ordering, labels)
  }
  values <- array(0, c(length(response), a$horizon + 1L, length(model$worlds)), dimnames = list(response, as.character(0:a$horizon), NULL))
  for (d in seq_along(model$worlds)) {
    w <- model$worlds[[d]]
    if (identification == "generalized") impact <- w$Sigma_u[, a$j] / sqrt(w$Sigma_u[a$j, a$j]) else {
      v <- t(chol(w$Sigma_u[ord, ord, drop = FALSE]))[, match(a$j, ord)]
      impact <- numeric(m); impact[ord] <- v
    }
    if (!is.null(shock_size)) {
      if (abs(impact[a$j]) < normalization_tol) stop("near-zero impact prevents requested normalization; vary normalization_tol for sensitivity")
      impact <- impact * shock_size / impact[a$j]
    }
    values[, , d] <- .recur_irf(w, impact, a$horizon)[a$r, , drop = FALSE]
  }
  .irf_result(values, model, impulse, response, identification, shock_size, list(ordering = ordering))
}
#' Compare two impulse response distributions
#' @param reference Reference JCERBGVAR_irf.
#' @param candidate Candidate JCERBGVAR_irf with identical selections and units.
#' @return Data frame with median RMSE and mean absolute interval endpoint difference.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 30)
#' f <- jcer_bgvar(s$data, s$W, draws = 10, seed = 1)
#' a <- jcer_irf(f, "US.r_short", "JP.y", horizon = 2)
#' jcer_compare_irf(a, a)
jcer_compare_irf <- function(reference, candidate) {
  for (key in c("response", "impulse", "horizon", "rate_unit", "identification", "shock_size", "ordering"))
    if (!identical(reference[[key]], candidate[[key]])) stop("IRF selections, normalization, identification and units must agree")
  delta <- reference$quantiles - candidate$quantiles
  data.frame(median_rmse = sqrt(mean(delta[3, , ]^2)), interval_endpoint_mae = mean(abs(delta[c(1, 5), , ])))
}
