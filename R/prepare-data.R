#' Prepare conditional country regressions
#' @param data Named list of numeric matrices with identical columns and sample indices.
#' @param W Named nonnegative country weight matrix, or a named list of such matrices per variable.
#' @param lags Positive lag order.
#' @param contemporaneous Include contemporaneous foreign variables.
#' @param normalize_weights Explicitly normalize weight rows after country selection.
#' @param rate_unit Interest rate input unit: percent points or decimal fractions.
#' @return A JCERBGVAR_data list containing country Y/Z matrices, extraction E/H maps, weights, names and units.
#' @export
#' @examples
#' s <- jcer_simulate_gvar(n_obs = 40)
#' d <- jcer_prepare_data(s$data, s$W)
jcer_prepare_data <- function(data, W, lags = 1, contemporaneous = TRUE,
                              normalize_weights = FALSE, rate_unit = c("percent", "decimal")) {
  lags <- .integer_arg(lags, "lags")
  rate_unit <- match.arg(rate_unit)
  if (!is.logical(contemporaneous) || length(contemporaneous) != 1L || is.na(contemporaneous) ||
      !is.logical(normalize_weights) || length(normalize_weights) != 1L || is.na(normalize_weights)) stop("logical options must be TRUE or FALSE")
  if (!is.list(data) || length(data) < 2 || !.named_unique(names(data))) stop("data needs at least two uniquely named countries")
  countries <- names(data)
  for (x in data) {
    if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) || !.named_unique(colnames(x)))
      stop("each country must be a finite numeric matrix with unique variable names; missing values are unsupported")
    if (!identical(dim(x), dim(data[[1]])) || !identical(colnames(x), colnames(data[[1]])))
      stop("v0.1 requires identical variable columns and sample sizes in every country")
    if (!identical(rownames(x), rownames(data[[1]])) || !identical(attr(x, "tsp"), attr(data[[1]], "tsp")))
      stop("country sample indices must agree")
  }
  vars <- colnames(data[[1]]); k <- length(vars); n <- length(data); m <- n * k
  if (any(grepl("\\.", c(countries, vars)))) stop("country and variable names cannot contain dots (reserved separator)")
  tt <- nrow(data[[1]])
  if (tt - lags < 5L) stop("series is too short after lags; need at least five usable observations")
  wl <- if (is.matrix(W)) setNames(rep(list(W), k), vars) else W
  if (!is.list(wl) || !.named_unique(names(wl)) || !setequal(names(wl), vars)) stop("variable-specific W must cover exactly the variables")
  wa <- lapply(wl[vars], .weights, countries = countries, normalize = normalize_weights)
  X <- do.call(cbind, data)
  labels <- unlist(lapply(countries, function(z) paste(z, vars, sep = ".")), use.names = FALSE)
  colnames(X) <- labels
  ii <- seq.int(lags + 1L, tt)
  blocks <- lapply(seq_len(n), function(i) {
    idx <- (i - 1L) * k + seq_len(k)
    E <- diag(m)[idx, , drop = FALSE]
    H <- matrix(0, k, m)
    for (v in seq_len(k)) H[v, seq.int(v, m, by = k)] <- wa[[v]]$W[i, ]
    foreign <- X %*% t(H)
    Z <- matrix(1, length(ii), 1L, dimnames = list(NULL, "intercept"))
    for (l in seq_len(lags)) Z <- cbind(Z, data[[i]][ii - l, , drop = FALSE])
    if (contemporaneous) Z <- cbind(Z, foreign[ii, , drop = FALSE])
    for (l in seq_len(lags)) Z <- cbind(Z, foreign[ii - l, , drop = FALSE])
    list(Y = data[[i]][ii, , drop = FALSE], Z = Z, E = E, H = H, index = idx)
  })
  names(blocks) <- countries
  structure(list(blocks = blocks, data = data, X = X, variables = vars, countries = countries,
                 labels = labels, weights = wa, lags = lags, contemporaneous = contemporaneous,
                 n_obs = length(ii), rate_unit = rate_unit, synthetic = isTRUE(attr(data, "synthetic"))), class = "JCERBGVAR_data")
}
