.haar <- function(n) {
  q <- qr(matrix(stats::rnorm(n * n), n, n))
  s <- sign(diag(qr.R(q))); s[s == 0] <- 1
  sweep(qr.Q(q), 2, s, "*")
}
#' Single-shock sign restrictions using Haar orthogonal rotations
#' @param model A JCERBGVAR fit.
#' @param restrictions Data frame with variable, horizon (integer >=0), sign (-1 or +1) columns.
#' @param impulse Shock normalization variable, constrained positive at impact.
#' @param response Response variable labels to retain.
#' @param horizon Last output horizon.
#' @param max_tries Maximum rotation attempts per posterior draw.
#' @param seed Seed for rotations.
#' @param shock_bp Rate impact in basis points; NULL keeps unit structural shock.
#' @param rotations_per_draw Accepted rotations required per posterior draw.
#' @param normalization_tol Minimum absolute impact before basis point scaling.
#' @return JCERBGVAR_irf. Errors if any posterior draw has insufficient accepted rotations; error condition contains attempts and accepted counts.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 40)
#' f <- jcer_bgvar(s$data, s$W, draws = 5, seed = 1)
#' r <- data.frame(variable = "US.r_short", horizon = 0L, sign = 1)
#' a <- jcer_irf_sign(f, r, "US.r_short", "JP.y", horizon = 3, seed = 2)
jcer_irf_sign <- function(model, restrictions, impulse, response, horizon = 20,
                           max_tries = 1000, seed = NULL, shock_bp = 25,
                           rotations_per_draw = 1, normalization_tol = 1e-8) {
  max_tries <- .integer_arg(max_tries, "max_tries"); nr <- .integer_arg(rotations_per_draw, "rotations_per_draw")
  shock_size <- NULL
  if (!is.null(shock_bp)) {
    if (length(shock_bp) != 1 || !is.finite(shock_bp) || shock_bp <= 0 || !grepl("\\.(r_short|r_long|stir|ltir)$", impulse)) stop("positive shock_bp requires a rate impulse")
    shock_size <- shock_bp / if (model$prepared$rate_unit == "percent") 100 else 10000
  }
  a <- .irf_args(model, impulse, response, horizon, shock_size, normalization_tol)
  labels <- model$prepared$labels
  if (!is.data.frame(restrictions) || !all(c("variable", "horizon", "sign") %in% names(restrictions)) ||
      anyNA(restrictions) || any(!restrictions$variable %in% labels) || !is.numeric(restrictions$horizon) ||
      any(!is.finite(restrictions$horizon)) || any(restrictions$horizon < 0 | restrictions$horizon != floor(restrictions$horizon)) ||
      any(!restrictions$sign %in% c(-1, 1))) stop("invalid sign restrictions: use variable, nonnegative integer horizon, sign +/-1")
  restrictions <- rbind(restrictions[, c("variable", "horizon", "sign")], data.frame(variable = impulse, horizon = 0L, sign = 1))
  hh <- max(c(a$horizon, restrictions$horizon)); ri <- match(restrictions$variable, labels)
  nd <- length(model$worlds); attempts <- accepted <- integer(nd)
  vals <- array(0, c(length(response), a$horizon + 1L, nd * nr), dimnames = list(response, as.character(0:a$horizon), NULL))
  .with_seed(seed, {
    for (d in seq_len(nd)) {
      w <- model$worlds[[d]]; P <- t(chol(w$Sigma_u))
      for (tr in seq_len(max_tries)) {
        attempts[d] <- tr
        impact <- drop(P %*% .haar(nrow(P))[, 1L])
        # Check unscaled signs; positive normalization cannot reverse restrictions.
        if (impact[a$j] <= normalization_tol) next
        path <- .recur_irf(w, impact, hh)
        if (all(path[cbind(ri, restrictions$horizon + 1L)] * restrictions$sign > 0)) {
          accepted[d] <- accepted[d] + 1L
          if (!is.null(shock_size)) path <- path * shock_size / impact[a$j]
          vals[, , (d - 1L) * nr + accepted[d]] <- path[a$r, seq_len(a$horizon + 1L), drop = FALSE]
          if (accepted[d] == nr) break
        }
      }
    }
  })
  if (any(accepted < nr)) {
    cond <- structure(list(message = paste("insufficient accepted rotations for", sum(accepted < nr), "posterior draws; no summary returned"),
                           call = NULL, attempts = attempts, accepted = accepted), class = c("JCERBGVAR_sign_failure", "error", "condition"))
    stop(cond)
  }
  .irf_result(vals, model, impulse, response, "single-shock sign restrictions", shock_size,
    list(restrictions = restrictions, attempts = attempts, accepted = accepted,
         acceptance_rate = sum(accepted) / sum(attempts), normalization_tol = normalization_tol,
         rotation_target = "country posterior draws kept equally weighted; conditional Haar rotations per draw; quantiles combine posterior and identification uncertainty"))
}
