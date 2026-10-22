#' Print a conditional Bayesian GVAR fit
#' @param x A fitted object.
#' @param ... Additional print arguments (unused).
#' @return x invisibly.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 30)
#' print(jcer_bgvar(s$data, s$W, draws = 5, seed = 1))
print.JCERBGVAR <- function(x, ...) {
  cat("preJCERBGVAR (development): conjugate MNIW independent country draws\n")
  cat(length(x$worlds), "usable /", x$requested_draws, "requested draws; rates:", x$prepared$rate_unit, "\n")
  cat(x$posterior_target, "\n")
  invisible(x)
}
#' Diagnose conditional posterior world draws
#' @param model A JCERBGVAR fit.
#' @return List with stability rate, minimum reciprocal condition number, rejection reasons, timings and posterior target.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 30)
#' f <- jcer_bgvar(s$data, s$W, draws = 5, seed = 1)
#' jcer_diagnostics(f)
jcer_diagnostics <- function(model) {
  if (!inherits(model, "JCERBGVAR")) stop("model must be a JCERBGVAR fit")
  list(method = model$method, approximation = model$approximation, burn_in = 0L,
       usable_draws = length(model$worlds), requested_draws = model$requested_draws,
       stable_rate = mean(vapply(model$worlds, function(w) w$stable, logical(1))),
       minimum_rcond = min(vapply(model$worlds, function(w) w$rcond, numeric(1))),
       rejected = model$rejected, rejection_reasons = model$rejection_reasons,
       posterior_target = model$posterior_target, rate_unit = model$prepared$rate_unit,
       timing = model$timing, object_bytes = as.numeric(object.size(model)),
       Rhat = NA_real_, ESS = length(model$worlds),
       MCMC_note = "Rhat not applicable to independent direct draws; ESS equals usable draws before sign rotations")
}
#' Summarize a conditional Bayesian GVAR fit
#' @param object A JCERBGVAR fit.
#' @param ... Unused.
#' @return Diagnostic list.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 30)
#' summary(jcer_bgvar(s$data, s$W, draws = 5, seed = 1))
summary.JCERBGVAR <- function(object, ...) jcer_diagnostics(object)
#' Plot selected impulse response distributions
#' @param x A JCERBGVAR_irf object.
#' @param ... Additional arguments passed to plot.
#' @return x invisibly. Draws median and 90 percent posterior/rotation quantiles.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 30)
#' f <- jcer_bgvar(s$data, s$W, draws = 10, seed = 1)
#' plot(jcer_irf(f, "US.r_short", "JP.y", horizon = 3))
plot.JCERBGVAR_irf <- function(x, ...) {
  old <- graphics::par(mfrow = c(length(x$response), 1), mar = c(3, 4, 3, 1))
  on.exit(graphics::par(old))
  for (r in seq_along(x$response)) {
    q <- matrix(x$quantiles[, r, ], 5)
    graphics::plot(x$horizon, q[3, ], type = "n", ylim = range(q), xlab = "Horizon", ylab = x$response[r],
      main = paste(if (x$synthetic) "SYNTHETIC DATA |" else "", x$identification), ...)
    graphics::polygon(c(x$horizon, rev(x$horizon)), c(q[1, ], rev(q[5, ])), col = "grey85", border = NA)
    graphics::lines(x$horizon, q[3, ], col = "navy", lwd = 2)
    graphics::abline(h = 0, lty = 3)
    graphics::mtext(paste(x$method, "| rates:", x$rate_unit, "| 90% quantiles"), side = 3, cex = 0.7)
  }
  invisible(x)
}
