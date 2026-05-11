#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "CellDynamicST v", utils::packageVersion("CellDynamicST"),
    " - Spatial Transcriptomic Analysis of Cell Dynamics\n",
    "Type ?CellDynamicST for an overview, or vignette('quickstart') to get started."
  )
}

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Set default options
  op <- options()
  op_cdst <- list(
    cdst.verbose = TRUE,
    cdst.parallel = FALSE,
    cdst.n_workers = 4L,
    cdst.seed = 42L
  )
  toset <- !(names(op_cdst) %in% names(op))
  if (any(toset)) options(op_cdst[toset])

  invisible()
}
