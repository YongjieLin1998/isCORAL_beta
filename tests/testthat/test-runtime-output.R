restore_test_omp <- function(value) {
  if (is.na(value)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = value)
}

test_that("OpenMP values distinguish unset, invalid and valid thread lists", {
  for (value in c("1", "8", " 2, 4 ", "0002")) {
    expect_identical(.isCORAL_normalize_omp_threads(value), value)
  }
  for (value in c(NA_character_, "", " ", "0", "-1", "1.5", "two", "2,0", "2,", "Inf")) {
    expect_identical(.isCORAL_normalize_omp_threads(value), NA_character_)
  }
  expect_identical(.isCORAL_normalize_omp_threads(NULL), NA_character_)
})

test_that("OpenMP restoration never leaves an empty environment variable", {
  original <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  on.exit(restore_test_omp(original), add = TRUE)
  for (value in c(NA_character_, "", " ", "0", "invalid")) {
    Sys.setenv(OMP_NUM_THREADS = "1")
    .isCORAL_restore_omp_threads(value)
    expect_identical(Sys.getenv("OMP_NUM_THREADS", unset = NA_character_), NA_character_)
  }
  for (value in c("4", "2,4")) {
    Sys.setenv(OMP_NUM_THREADS = "1")
    .isCORAL_restore_omp_threads(value)
    expect_identical(Sys.getenv("OMP_NUM_THREADS"), value)
  }
})

test_that("CORALembedding restores OpenMP state even after an early error", {
  original <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  on.exit(restore_test_omp(original), add = TRUE)
  for (value in c(NA_character_, "", "invalid", "4", "2,4")) {
    restore_test_omp(value)
    # Invalid object fails after the temporary thread environment is installed.
    expect_error(CORALembedding(NULL))
    expect_identical(
      Sys.getenv("OMP_NUM_THREADS", unset = NA_character_),
      .isCORAL_normalize_omp_threads(value)
    )
  }
})

test_that("scan child processes receive thread limits before package loading", {
  skip_if_not_installed("callr")
  original <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  on.exit(restore_test_omp(original), add = TRUE)
  Sys.setenv(OMP_NUM_THREADS = "invalid")
  child_env <- .isCORAL_callr_thread_env()
  thread_names <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS")
  expect_true(all(child_env[thread_names] == "1"))
  safe_env <- callr::rcmd_safe_env()
  other_names <- setdiff(names(safe_env), thread_names)
  expect_identical(child_env[other_names], safe_env[other_names])
  observed <- callr::r(
    function(keys) Sys.getenv(keys), args = list(keys = thread_names),
    env = child_env, user_profile = FALSE
  )
  expect_identical(unname(observed), rep("1", length(thread_names)))
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "invalid")
})

test_that("quiet Seurat output preserves results, warnings, errors and sinks", {
  result <- NULL
  expect_output(expect_message({
    result <- .isCORAL_with_seurat_output({
      cat("routine stdout\n")
      message("routine message")
      list(value = 42L)
    })
  }, NA), NA)
  expect_identical(result, list(value = 42L))
  expect_warning(.isCORAL_with_seurat_output(warning("important warning")), "important warning")
  sinks_before <- sink.number()
  expect_error(.isCORAL_with_seurat_output({
    cat("partial stdout\n")
    stop("important failure")
  }), "important failure")
  expect_identical(sink.number(), sinks_before)
  expect_output(.isCORAL_with_seurat_output(cat("debug stdout"), verbose = TRUE), "debug stdout")
  expect_message(.isCORAL_with_seurat_output(message("debug message"), verbose = TRUE), "debug message")
})

test_that("quiet object updates retain the same Seurat data", {
  skip_if_not_installed("SeuratObject")
  counts <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  object <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  quiet <- .isCORAL_with_seurat_output(SeuratObject::UpdateSeuratObject(object))
  baseline <- suppressMessages(SeuratObject::UpdateSeuratObject(object))
  expect_identical(quiet, baseline)
})

test_that("public analysis entries default to quiet routine output", {
  expect_identical(formals(isCORAL)$verbose, FALSE)
  expect_identical(formals(CORALembedding)$verbose, FALSE)
  expect_identical(formals(FindMetaclone)$verbose, FALSE)
  expect_identical(formals(FindMetaclones), formals(FindMetaclone))
})

test_that("quiet and verbose Metaclone runs preserve labels and graphs", {
  skip_if_not_installed("Seurat")
  set.seed(91)
  counts <- matrix(stats::rpois(400, 4), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("c", 1:40)))
  object <- Seurat::CreateSeuratObject(counts = counts)
  coordinates <- matrix(stats::rnorm(120), ncol = 3,
                        dimnames = list(colnames(object), paste0("CORALclone_", 1:3)))
  object[["CORAL_clone_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = coordinates, key = "CORALclone_", assay = "RNA"
  )
  quiet <- NULL
  expect_output(expect_message({
    quiet <- FindMetaclone(object, dims = 1:3, algorithm = 1, random_seed = 29)
  }, NA), NA)
  loud <- NULL
  invisible(capture.output(loud <- suppressMessages(
    FindMetaclone(object, dims = 1:3, algorithm = 1, random_seed = 29, verbose = TRUE)
  )))
  expect_identical(quiet$metaclone, loud$metaclone)
  expect_equal(quiet@graphs, loud@graphs)
})
