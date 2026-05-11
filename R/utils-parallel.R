#' @title Parallel Execution Utilities
#' @description Internal wrappers for parallel execution using future or
#'   BiocParallel backends.
#' @name parallel
#' @keywords internal
NULL

#' Apply a function in parallel over a list
#'
#' Provides a unified interface for parallel execution. When parallel = TRUE,
#' uses future.apply if available, otherwise falls back to sequential lapply.
#'
#' @param x A list or vector to iterate over.
#' @param fun A function to apply to each element.
#' @param parallel Logical. Use parallel execution. Default: FALSE.
#' @param n_workers Integer. Number of parallel workers. Default: 4.
#' @param ... Additional arguments passed to fun.
#'
#' @return A list of results.
#' @keywords internal
cdst_lapply <- function(x, fun, parallel = FALSE, n_workers = 4L, ...) {
  if (!parallel) {
    return(lapply(x, fun, ...))
  }

  if (requireNamespace("future.apply", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_workers)
    future.apply::future_lapply(x, fun, ..., future.seed = TRUE)
  } else {
    cli::cli_warn(c(
      "Parallel execution requested but {.pkg future.apply} not installed.",
      "i" = "Falling back to sequential execution.",
      "i" = "Install with: {.code install.packages(\"future.apply\")}"
    ))
    lapply(x, fun, ...)
  }
}


#' Set Up Parallel Backend
#'
#' Configures the parallel execution backend for CellDynamicST. Uses the
#' \pkg{future} package to set up multisession parallelism. If \pkg{future}
#' is not installed, a warning is issued and sequential execution is used.
#'
#' @param n_cores Integer. Number of cores to use. Default: 1 (sequential).
#' @param verbose Logical. Print informational messages. Default: TRUE.
#'
#' @return Invisibly returns the previous plan, which can be used to restore
#'   the original state.
#'
#' @examples
#' \dontrun{
#' old_plan <- setup_parallel(n_cores = 4)
#' # ... run parallel operations ...
#' future::plan(old_plan)  # restore
#' }
#'
#' @export
setup_parallel <- function(n_cores = 1L, verbose = TRUE) {
  checkmate::assert_int(n_cores, lower = 1L)

  if (n_cores == 1L) {
    if (verbose) cli::cli_inform("Using sequential execution (n_cores = 1)")
    return(invisible(NULL))
  }

  if (!requireNamespace("future", quietly = TRUE)) {
    cli::cli_warn(c(
      "{.pkg future} not installed. Using sequential execution.",
      "i" = "Install with: {.code install.packages(\"future\")}"
    ))
    return(invisible(NULL))
  }

  old_plan <- future::plan()
  future::plan(future::multisession, workers = n_cores)

  if (verbose) {
    cli::cli_inform("Parallel backend: {n_cores} workers (future::multisession)")
  }

  invisible(old_plan)
}
