.prior_mniw <- function(z, k, prior_control) {
  defaults <- list(B0 = matrix(0, z, k), V0 = diag(1, z), S0 = diag(0.1, k), nu0 = k + 2)
  if (!is.list(prior_control) || any(!names(prior_control) %in% names(defaults))) stop("unknown prior_control item")
  pr <- utils::modifyList(defaults, prior_control)
  if (!is.matrix(pr$B0) || !identical(dim(pr$B0), c(z, k)) || any(!is.finite(pr$B0))) stop("B0 has wrong dimensions or nonfinite values")
  if (!identical(dim(pr$V0), c(z, z)) || !identical(dim(pr$S0), c(k, k))) stop("prior covariance dimensions do not match")
  .spd(pr$V0, "V0"); .spd(pr$S0, "S0")
  if (length(pr$nu0) != 1 || !is.finite(pr$nu0) || pr$nu0 <= k - 1) stop("nu0 must exceed k - 1")
  pr
}
.mniw_posterior <- function(Y, Z, pr) {
  Q0 <- .spd_solve(pr$V0, diag(nrow(pr$V0)))
  Qn <- Q0 + crossprod(Z)
  Bn <- .spd_solve(Qn, Q0 %*% pr$B0 + crossprod(Z, Y))
  D <- Bn - pr$B0; residual <- Y - Z %*% Bn
  Sn <- pr$S0 + crossprod(residual) + crossprod(D, Q0 %*% D)
  list(Bn = Bn, Vn = .spd_solve(Qn, diag(nrow(Qn))), Sn = (Sn + t(Sn)) / 2,
       nun = pr$nu0 + nrow(Y), Qn = Qn, prior = pr)
}
.riw <- function(S, nu, wishart_scale = .spd_solve(S, diag(nrow(S)))) {
  precision <- stats::rWishart(1L, nu, wishart_scale)[, , 1L, drop = FALSE]
  precision <- matrix(precision, nrow(S), nrow(S))
  .spd_solve(precision, diag(nrow(S)))
}
.country_draws <- function(task) {
  .with_seed(task$seed, {
    post <- .mniw_posterior(task$Y, task$Z, task$prior)
    k <- ncol(task$Y); z <- ncol(task$Z)
    Bs <- array(0, c(z, k, task$draws)); Ss <- array(0, c(k, k, task$draws))
    Lv <- t(chol(post$Vn))
    wishart_scale <- .spd_solve(post$Sn, diag(k))
    for (d in seq_len(task$draws)) {
      S <- .riw(post$Sn, post$nun, wishart_scale)
      Bs[, , d] <- post$Bn + Lv %*% matrix(stats::rnorm(z * k), z, k) %*% chol(S)
      Ss[, , d] <- S
    }
    list(B = Bs, Sigma = Ss, posterior = post)
  })
}
#' Fit a conditional Bayesian GVAR with independent conjugate draws
#' @param data Named country matrices.
#' @param W Country weights (matrix or named variable-specific list).
#' @param lags Positive lag order.
#' @param method Only 'conjugate' is validated in this release.
#' @param prior Only 'mniw' is supported.
#' @param draws Number of independent posterior draws; no burn-in.
#' @param seed Optional reproducible seed; local RNG state is restored.
#' @param cores Country worker count; reproducible across worker counts.
#' @param contemporaneous Include contemporaneous foreign regressors.
#' @param normalize_weights Explicit weight normalization.
#' @param rate_unit Interest rate unit, 'percent' or 'decimal'.
#' @param prior_control Optional B0, V0, S0, nu0, common to all countries.
#' @param stability Either 'keep' all draws or 'condition' on a stable world system.
#' @return JCERBGVAR object with prepared data, country draws, world systems, timings, assumptions and rejection counts.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 40)
#' f <- jcer_bgvar(s$data, s$W, draws = 20, seed = 1)
#' summary(f)
jcer_bgvar <- function(data, W, lags = 1, method = "conjugate", prior = "mniw",
                        draws = 2000, seed = NULL, cores = 1, contemporaneous = TRUE,
                        normalize_weights = FALSE, rate_unit = c("percent", "decimal"),
                        prior_control = list(), stability = c("keep", "condition")) {
  if (!identical(method, "conjugate") || !identical(prior, "mniw")) stop("only method='conjugate', prior='mniw' is validated; VB/hybrid are research work")
  draws <- .integer_arg(draws, "draws"); cores <- .integer_arg(cores, "cores")
  stability <- match.arg(stability); rate_unit <- match.arg(rate_unit)
  start <- proc.time()[[3L]]
  pre <- jcer_prepare_data(data, W, lags, contemporaneous, normalize_weights, rate_unit)
  tpre <- proc.time()[[3L]] - start
  tasks <- .with_seed(seed, lapply(seq_along(pre$blocks), function(i) {
    block <- pre$blocks[[i]]
    list(Y = block$Y, Z = block$Z, prior = .prior_mniw(ncol(block$Z), ncol(block$Y), prior_control),
         seed = sample.int(.Machine$integer.max, 1), draws = draws)
  }))
  t0 <- proc.time()[[3L]]
  if (cores > 1L) {
    cl <- parallel::makePSOCKcluster(min(cores, length(tasks)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
    countries <- parallel::parLapply(cl, tasks, .country_draws)
  } else countries <- lapply(tasks, .country_draws)
  names(countries) <- pre$countries
  tcountry <- proc.time()[[3L]] - t0; t0 <- proc.time()[[3L]]
  worlds <- vector("list", draws); invalid <- character(draws); keep <- rep(TRUE, draws)
  for (d in seq_len(draws)) {
    w <- tryCatch(.world_system(lapply(countries, function(x) matrix(x$B[, , d], nrow(x$B), dim(x$B)[2])),
                               lapply(countries, function(x) matrix(x$Sigma[, , d], dim(x$Sigma)[1], dim(x$Sigma)[2])), pre), error = identity)
    if (inherits(w, "error")) { keep[d] <- FALSE; invalid[d] <- conditionMessage(w) } else {
      worlds[[d]] <- w
      if (stability == "condition" && !w$stable) { keep[d] <- FALSE; invalid[d] <- "unstable world" }
    }
  }
  if (!all(keep)) warning(sum(!keep), " / ", draws, " draws excluded; inspect diagnostics and changed posterior target")
  if (!any(keep)) stop("no usable world draws")
  target <- "product of country conditional MNIW posteriors; weak foreign exogeneity assumed; not full simultaneous-world posterior"
  if (stability == "condition") target <- paste(target, "conditioned on stable world dynamics", sep = "; ")
  structure(list(prepared = pre, countries = countries, worlds = worlds[keep], draw_ids = which(keep),
    method = method, prior = prior, seed = seed, cores = cores, requested_draws = draws,
    rejection_reasons = table(invalid[nzchar(invalid)]), rejected = sum(!keep), stability = stability,
    posterior_target = target, approximation = "exact independent draws for each conditional country posterior",
    timing = c(preprocessing = tpre, country_estimation = tcountry, world_construction = proc.time()[[3L]] - t0,
               total = proc.time()[[3L]] - start)), class = "JCERBGVAR")
}
