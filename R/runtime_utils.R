# Runtime helpers for thread environments and routine console output.

.isCORAL_normalize_omp_threads <- function(value) {
  # OpenMP accepts positive integers, including lists for nested parallelism.
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[[:space:]]*[0-9]+[[:space:]]*(,[[:space:]]*[0-9]+[[:space:]]*)*$", value)) {
    return(NA_character_)
  }
  counts <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (any(!is.finite(counts)) || any(counts < 1 | counts > .Machine$integer.max)) {
    return(NA_character_)
  }
  value
}

.isCORAL_restore_omp_threads <- function(value) {
  value <- .isCORAL_normalize_omp_threads(value)
  if (is.na(value)) {
    # An unset variable is not equivalent to OMP_NUM_THREADS="".
    Sys.unsetenv("OMP_NUM_THREADS")
  } else {
    Sys.setenv(OMP_NUM_THREADS = value)
  }
  invisible(value)
}

.isCORAL_callr_thread_env <- function() {
  env <- callr::rcmd_safe_env()
  env[c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
        "RCPP_PARALLEL_NUM_THREADS")] <- "1"
  env
}

.isCORAL_with_seurat_output <- function(expr, verbose = FALSE) {
  if (isTRUE(verbose)) return(force(expr))
  result <- NULL
  # Suppress routine messages/stdout only; warnings and errors still propagate.
  invisible(utils::capture.output(result <- suppressMessages(force(expr))))
  result
}
