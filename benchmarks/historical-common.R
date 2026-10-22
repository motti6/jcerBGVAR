# Reproducible data preparation and public BGVAR API helpers; no data bundled.
historical_data <- function(n_countries = 4L, include_long_rate = FALSE,
                            input_rate_unit = "decimal", transformation = "growth") {
  if (!requireNamespace("BGVAR", quietly = TRUE)) stop("Install optional BGVAR dependency")
  e <- new.env(parent = emptyenv()); utils::data("eerData", package = "BGVAR", envir = e)
  vars <- c("y", "Dp", "stir", if (include_long_rate) "ltir")
  candidates <- c("US", "JP", "EA", "KR", setdiff(names(e$eerData), c("US", "JP", "EA", "KR")))
  eligible <- candidates[vapply(e$eerData[candidates], function(x) all(vars %in% colnames(x)) && all(is.finite(x[, vars])), logical(1))]
  if (n_countries > length(eligible)) stop("not enough countries with selected variables")
  cc <- eligible[seq_len(n_countries)]
  if (!input_rate_unit %in% c("decimal", "percent")) stop("explicit rate unit required")
  factor <- if (input_rate_unit == "decimal") 100 else 1
  d <- lapply(e$eerData[cc], function(x) {
    idx <- if (inherits(x, "zoo")) as.character(zoo::index(x)) else if (is.ts(x)) as.character(zoo::as.yearqtr(stats::time(x))) else rownames(x)
    x <- matrix(as.numeric(x), nrow(x), ncol(x), dimnames = dimnames(x))
    if (is.null(idx)) idx <- paste0(1995L + (seq_len(nrow(x)) - 1L) %/% 4L, "Q", (seq_len(nrow(x)) - 1L) %% 4L + 1L)
    if (transformation == "growth") {
      z <- cbind(y = 100 * diff(x[, "y"]), Dp = 100 * x[-1, "Dp"], stir = factor * x[-1, "stir"])
      if (include_long_rate) z <- cbind(z, ltir = factor * x[-1, "ltir"])
      rownames(z) <- idx[-1]
    } else if (transformation == "levels") {
      z <- x[, vars, drop = FALSE]; z[, "Dp"] <- z[, "Dp"] * 100
      z[, intersect(c("stir", "ltir"), vars)] <- z[, intersect(c("stir", "ltir"), vars), drop = FALSE] * factor
      rownames(z) <- idx
    } else stop("unknown transformation")
    colnames(z) <- vars
    z
  })
  W <- e$W.trade0012[cc, cc, drop = FALSE]; original <- rowSums(W); W <- W / original
  list(data = d, W = W, countries = cc, variables = vars, retained_trade_mass = original,
       source_package = as.character(utils::packageVersion("BGVAR")), rate_unit = "percent",
       input_rate_unit = input_rate_unit, transformation = transformation,
       period = range(rownames(d[[1]])), excluded = setdiff(names(e$eerData), cc))
}

historical_bgvar_fit <- function(d, draws, burnin, seed, prior = "NG", SV = FALSE, lags = 1) {
  set.seed(seed)
  BGVAR::bgvar(Data = d$data, W = d$W, plag = lags, draws = draws, burnin = burnin,
    prior = prior, SV = SV, thin = 1, trend = FALSE, eigen = FALSE,
    expert = list(applyfun = lapply, cores = NULL), verbose = FALSE)
}
historical_bgvar_irf <- function(fit, horizon = 12, shock_size = .25) {
  shock <- BGVAR::get_shockinfo("girf"); shock$shock <- "US.stir"; shock$scale <- shock_size
  BGVAR::irf(fit, n.ahead = horizon, shockinfo = shock, quantiles = c(.05, .16, .5, .84, .95),
             expert = list(save.store = TRUE, cores = NULL), verbose = FALSE)
}
historical_irf_array <- function(irf) {
  arr <- irf$IRF_store
  out <- array(arr[, 1L, , ], c(dim(arr)[1], dim(arr)[3], dim(arr)[4]))
  dimnames(out) <- list(dimnames(arr)[[1]], dimnames(arr)[[3]], NULL)
  out
}
stable_bgvar_subset <- function(fit) {
  idx <- which(as.numeric(fit$stacked.results$F.eigen) < 1)
  if (!length(idx)) stop("no strictly stable BGVAR draws")
  out <- fit
  for (key in c("S_large", "Ginv_large", "A_large"))
    out$stacked.results[[key]] <- fit$stacked.results[[key]][, , idx, drop = FALSE]
  out$stacked.results$F_large <- fit$stacked.results$F_large[, , , idx, drop = FALSE]
  out$stacked.results$F.eigen <- as.numeric(fit$stacked.results$F.eigen)[idx]
  out$args$thindraws <- length(idx)
  out
}
