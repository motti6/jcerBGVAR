.integer_arg <- function(x, name, lower = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) || x < lower)
    stop(name, " must be an integer >= ", lower, call. = FALSE)
  as.integer(x)
}
.named_unique <- function(x) !is.null(x) && all(nzchar(x)) && !anyNA(x) && !anyDuplicated(x)
.with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  seed <- .integer_arg(seed, "seed", 0L)
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  force(code)
}
.spd <- function(x, name = "matrix") {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x) || any(!is.finite(x)) ||
      !isTRUE(all.equal(x, t(x), tolerance = 1e-10, check.attributes = FALSE)))
    stop(name, " must be symmetric positive definite", call. = FALSE)
  tryCatch(chol(x), error = function(e) stop(name, " must be positive definite", call. = FALSE))
}
.spd_solve <- function(x, b) {
  u <- .spd(x)
  backsolve(u, forwardsolve(t(u), b))
}
.weights <- function(W, countries, normalize) {
  if (!is.matrix(W) || !is.numeric(W) || nrow(W) != ncol(W) ||
      !.named_unique(rownames(W)) || !.named_unique(colnames(W)) ||
      !setequal(rownames(W), colnames(W)) || !all(countries %in% rownames(W)) ||
      any(!is.finite(W)) || any(W < 0)) stop("W must be a finite nonnegative named square matrix covering all countries")
  if (any(abs(diag(W[rownames(W), rownames(W), drop = FALSE])) > 1e-12)) stop("W self weights must be zero")
  excluded <- setdiff(rownames(W), countries)
  W <- W[countries, countries, drop = FALSE]
  s <- rowSums(W)
  if (any(s <= 0)) stop("W has a zero foreign-weight row after country selection")
  if (!normalize && any(abs(s - 1) > 1e-8)) stop("W row sums must be one; set normalize_weights = TRUE explicitly")
  if (normalize) W <- W / s
  list(W = W, excluded = excluded, original_row_sums = s, renormalized = normalize && any(abs(s - 1) > 1e-8))
}
