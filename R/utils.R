# Internal utilities ---------------------------------------------------------

.check_scalar <- function(x, name, lower = -Inf, upper = Inf, strict = FALSE,
                          integer = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) {
    ok <- if (strict) x > lower else x >= lower
    ok <- ok && x <= upper
  }
  if (ok && integer) ok <- abs(x - round(x)) < sqrt(.Machine$double.eps)
  if (!ok) {
    rng <- if (strict) sprintf("> %s", lower) else sprintf(">= %s", lower)
    if (is.finite(upper)) rng <- paste(rng, sprintf("and <= %s", upper))
    stop(sprintf("`%s` must be a single finite %s %s.", name,
                 if (integer) "integer" else "number", rng), call. = FALSE)
  }
  invisible(x)
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(x)
}

.as_numeric_matrix <- function(x, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x) && is.numeric(x)) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix or data frame.", name),
         call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` contains missing or non-finite values.", name),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

# Apply FUN over X, optionally in parallel through future.apply.
# FUN must not use random numbers (all randomness is drawn beforehand so
# results do not depend on the parallel backend).
.cpf_lapply <- function(X, FUN, parallel = FALSE, progress = FALSE,
                        label = NULL) {
  n <- length(X)
  if (n == 0L) return(list())
  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("`parallel = TRUE` requires the 'future.apply' package. Install it ",
           "and choose a backend with future::plan(), e.g. ",
           "future::plan(future::multisession).", call. = FALSE)
    }
    if (isTRUE(progress) && !is.null(label)) {
      message(label, " (", n, " fits, in parallel)...")
    }
    return(future.apply::future_lapply(X, FUN, future.seed = NULL))
  }
  if (isTRUE(progress) && !is.null(label)) message(label, " (", n, " fits)...")
  pb <- if (isTRUE(progress) && interactive() && n > 1L) {
    utils::txtProgressBar(min = 0, max = n, style = 3)
  } else {
    NULL
  }
  out <- vector("list", n)
  for (i in seq_len(n)) {
    out[[i]] <- FUN(X[[i]])
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
  }
  if (!is.null(pb)) close(pb)
  out
}

# Run `expr` with a temporary RNG seed, restoring the caller's RNG state.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  .check_scalar(seed, "seed")
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expr)
}

.fmt <- function(x, digits = 4) formatC(x, digits = digits, format = "f")
