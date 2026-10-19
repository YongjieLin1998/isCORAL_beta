make_test_iscoral_object <- function(n_cells = 30L, n_genes = 12L) {
  set.seed(101)
  counts <- matrix(
    stats::rpois(n_cells * n_genes, lambda = 4),
    nrow = n_genes,
    dimnames = list(
      paste0("gene_", seq_len(n_genes)),
      paste0("cell_", seq_len(n_cells))
    )
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  Seurat::NormalizeData(object, verbose = FALSE)
}

test_that("historical Metaclone formal parameters remain compatible", {
  expect_identical(formals(FindMetaclone), formals(FindMetaclones))
  expect_identical(formals(CORALembedding)$max_cells_lasso, NULL)
  expect_identical(formals(CORALembedding)$random_seed, NULL)
})

test_that("Leiden dependency policy resolves explicitly to Leiden or Louvain", {
  resolver <- getFromNamespace(".resolve_metaclone_algorithm", "isCORAL")

  leiden_result <- resolver(
    algorithm = 4,
    verbose = FALSE,
    leidenbase_available = TRUE,
    seurat_version = "5.2.0"
  )
  expect_identical(leiden_result$algorithm_requested, 4)
  expect_identical(leiden_result$algorithm_effective, 4)
  expect_identical(leiden_result$clustering_engine, "Leiden (leidenbase)")

  fallback_result <- expect_warning(
    resolver(
      algorithm = 4,
      verbose = FALSE,
      leidenbase_available = FALSE,
      seurat_version = "5.2.0"
    ),
    "does not install optional packages at runtime"
  )
  expect_identical(fallback_result$algorithm_requested, 4)
  expect_identical(fallback_result$algorithm_effective, 1)
  expect_identical(
    fallback_result$clustering_engine,
    "Louvain fallback (leidenbase unavailable)"
  )

  legacy_seurat_result <- expect_warning(
    resolver(
      algorithm = 4,
      verbose = FALSE,
      leidenbase_available = TRUE,
      seurat_version = "5.1.0"
    ),
    "Seurat >= 5.2.0"
  )
  expect_identical(legacy_seurat_result$algorithm_requested, 4)
  expect_identical(legacy_seurat_result$algorithm_effective, 1)
  expect_identical(
    legacy_seurat_result$clustering_engine,
    "Louvain fallback (Seurat < 5.2.0)"
  )
})

test_that("leidenbase is the declared optional Leiden dependency", {
  package_description <- utils::packageDescription("isCORAL")
  suggests <- trimws(strsplit(package_description$Suggests, ",", fixed = TRUE)[[1]])

  expect_true("leidenbase" %in% suggests)
  expect_false("leiden" %in% suggests)
})

test_that("FindMetaclone writes aligned metadata and a run record", {
  skip_if_not_installed("Seurat")
  object <- make_test_iscoral_object()
  set.seed(102)
  clone_coordinates <- matrix(
    stats::rnorm(ncol(object) * 3L),
    nrow = ncol(object),
    dimnames = list(colnames(object), paste0("CORALclone_", 1:3))
  )
  object[["CORAL_clone_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = clone_coordinates,
    key = "CORALclone_",
    assay = Seurat::DefaultAssay(object)
  )

  result <- expect_warning(
    FindMetaclone(
      object,
      resolution = 0.2,
      graph_prefix = "compat",
      cluster_col_name = "metaclone_test",
      algorithm = 1,
      random_seed = 17,
      verbose = FALSE
    ),
    "historical default"
  )

  expect_true("metaclone_test" %in% colnames(result[[]]))
  expect_identical(
    as.character(Seurat::Idents(result)),
    as.character(result$metaclone_test)
  )
  expect_true(all(c("compat_nn", "compat_snn") %in% names(result@graphs)))
  expect_identical(
    result@misc$isCORAL_metaclone$latest$effective_parameters$algorithm,
    1
  )
  expect_identical(
    result@misc$isCORAL_metaclone$latest$requested_parameters$dims,
    1:30
  )
  expect_identical(
    result@misc$isCORAL_metaclone$latest$effective_parameters$dims,
    1:3
  )
  expect_identical(
    result@misc$isCORAL_metaclone$latest$outputs$cluster_metadata_column,
    "metaclone_test"
  )
})

test_that("historical FindMetaclones entry point forwards to FindMetaclone", {
  skip_if_not_installed("Seurat")
  object <- make_test_iscoral_object()
  set.seed(103)
  clone_coordinates <- matrix(
    stats::rnorm(ncol(object) * 3L),
    nrow = ncol(object),
    dimnames = list(colnames(object), paste0("CORALclone_", 1:3))
  )
  object[["CORAL_clone_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = clone_coordinates,
    key = "CORALclone_",
    assay = Seurat::DefaultAssay(object)
  )

  canonical <- FindMetaclone(
    object, 0.2, "CORAL_clone_embed", 1:3, "legacy", "metaclone",
    1, 29, FALSE
  )
  legacy <- expect_warning(
    FindMetaclones(
      object, 0.2, "CORAL_clone_embed", 1:3, "legacy", "metaclone",
      1, 29, FALSE
    ),
    "deprecated"
  )

  expect_identical(
    as.character(canonical$metaclone),
    as.character(legacy$metaclone)
  )
  expect_identical(
    legacy@misc$isCORAL_metaclone$latest$entry_point,
    "FindMetaclones"
  )
  expect_match(
    legacy@misc$isCORAL_metaclone$latest$call_text,
    "^FindMetaclones"
  )
})

test_that("CORALembedding retains legacy cache fields and records parameters", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("irlba")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("foreach")
  skip_if_not_installed("matrixStats")

  object <- make_test_iscoral_object(n_cells = 20L, n_genes = 10L)
  cells <- colnames(object)
  probability <- Matrix::sparseMatrix(
    i = seq_along(cells),
    j = seq_along(cells),
    x = 1,
    dims = c(length(cells), length(cells)),
    dimnames = list(cells, NULL)
  )
  gene_statistics <- data.frame(
    omega_square_infer = seq(1, 0.1, length.out = nrow(object)),
    row.names = rownames(object)
  )
  object@misc$CORAL_scan_results <- list(
    binsize_200 = list(
      S_stab = probability,
      ss_filter_average = gene_statistics
    )
  )

  foreach_globals <- get(
    ".foreachGlobals",
    envir = asNamespace("foreach"),
    inherits = FALSE
  )
  initial_backend <- NULL
  if (foreach::getDoParRegistered()) {
    initial_backend <- list(
      fun = get("fun", envir = foreach_globals, inherits = FALSE),
      data = get("data", envir = foreach_globals, inherits = FALSE),
      info = get("info", envir = foreach_globals, inherits = FALSE)
    )
  }
  restore_initial_backend <- function() {
    if (is.null(initial_backend)) {
      current_bindings <- intersect(
        c("fun", "data", "info"),
        ls(envir = foreach_globals, all.names = TRUE)
      )
      if (length(current_bindings) > 0L) {
        rm(list = current_bindings, envir = foreach_globals)
      }
    } else {
      foreach::setDoPar(
        initial_backend$fun,
        initial_backend$data,
        initial_backend$info
      )
    }
  }
  on.exit(restore_initial_backend(), add = TRUE)

  backend_bindings <- intersect(
    c("fun", "data", "info"),
    ls(envir = foreach_globals, all.names = TRUE)
  )
  if (length(backend_bindings) > 0L) {
    rm(list = backend_bindings, envir = foreach_globals)
  }
  expect_false(foreach::getDoParRegistered())
  backend_before <- list(
    name = foreach::getDoParName(),
    version = foreach::getDoParVersion(),
    workers = foreach::getDoParWorkers()
  )
  result <- CORALembedding(
    object,
    binsize_key = "binsize_200",
    n_dims = 2,
    top_n_genes = 8,
    max_cells_lasso = NULL,
    alpha_elasticnet = 0.05,
    workers_glmnet = 1,
    random_seed = 41
  )
  cache <- result@misc$isCORAL_embedding

  legacy_fields <- c(
    "P_norm", "original_MDS_coords", "Y_ref_raw", "valid_cells",
    "top_genes", "source_binsize", "lasso_model"
  )
  expect_true(all(legacy_fields %in% names(cache)))
  expect_identical(cache$iscoral_version, "1.0.0")
  expect_identical(cache$release_status, "beta")
  expect_identical(cache$run_parameters$n_dims, 2L)
  expect_identical(cache$run_parameters$top_n_genes, 8L)
  expect_null(cache$run_parameters$max_cells_lasso)
  expect_identical(cache$run_parameters$random_seed, 41L)
  expect_false(cache$effective_parameters$subsampling_applied)
  expect_identical(cache$effective_parameters$n_valid_cells, 20L)
  expect_identical(cache$effective_parameters$n_training_cells, 20L)
  expect_identical(cache$effective_parameters$effective_dimensions, 2L)
  expect_identical(cache$effective_parameters$glmnet$nfolds_effective, 5L)
  expect_identical(
    cache$provenance$probability_dimname_action,
    "copied_from_rows"
  )
  expect_identical(length(cache$reproducibility$training_cells), 20L)
  expect_identical(length(cache$reproducibility$cv_foldid), 20L)
  expect_identical(length(unique(cache$reproducibility$cv_foldid)), 5L)
  expect_identical(
    list(
      name = foreach::getDoParName(),
      version = foreach::getDoParVersion(),
      workers = foreach::getDoParWorkers()
    ),
    backend_before
  )
  expect_true(all(c("CORAL_cell_embed", "CORAL_clone_embed") %in%
                    names(result@reductions)))

  foreach::registerDoSEQ()
  sequential_backend <- list(
    name = foreach::getDoParName(),
    version = foreach::getDoParVersion(),
    workers = foreach::getDoParWorkers()
  )
  repeated_result <- CORALembedding(
    object,
    binsize_key = "binsize_200",
    n_dims = 2,
    top_n_genes = 8,
    max_cells_lasso = NULL,
    alpha_elasticnet = 0.05,
    workers_glmnet = 1,
    random_seed = 41
  )
  repeated_cache <- repeated_result@misc$isCORAL_embedding
  expect_identical(
    list(
      name = foreach::getDoParName(),
      version = foreach::getDoParVersion(),
      workers = foreach::getDoParWorkers()
    ),
    sequential_backend
  )
  expect_identical(
    repeated_cache$reproducibility$cv_foldid,
    cache$reproducibility$cv_foldid
  )
  expect_equal(
    repeated_cache$Y_ref_raw,
    cache$Y_ref_raw,
    tolerance = 1e-10
  )

  cache_file <- tempfile(fileext = ".rds")
  saveRDS(cache, cache_file)
  restored_cache <- readRDS(cache_file)
  expect_identical(restored_cache$reproducibility$cv_foldid,
                   cache$reproducibility$cv_foldid)
})

test_that("CORALembedding records finite training-cell subsampling", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("irlba")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("matrixStats")

  object <- make_test_iscoral_object(n_cells = 20L, n_genes = 10L)
  cells <- colnames(object)
  probability <- Matrix::Diagonal(n = length(cells), x = 1)
  dimnames(probability) <- list(cells, cells)
  gene_statistics <- data.frame(
    omega_square_infer = seq(1, 0.1, length.out = nrow(object)),
    row.names = rownames(object)
  )
  object@misc$coral_scan_results <- list(
    binsize_200 = list(
      S_stab = probability,
      ss_filter_average = gene_statistics
    )
  )

  result <- CORALembedding(
    object,
    n_dims = 2,
    top_n_genes = 8,
    max_cells_lasso = 15,
    workers_glmnet = 1,
    random_seed = 43
  )
  cache <- result@misc$isCORAL_embedding

  expect_identical(cache$effective_parameters$scan_misc_key, "coral_scan_results")
  expect_true(cache$effective_parameters$subsampling_applied)
  expect_identical(cache$effective_parameters$n_training_cells, 15L)
  expect_identical(length(cache$reproducibility$training_cells), 15L)
  expect_identical(length(cache$reproducibility$cv_foldid), 15L)
  seurat_object_v5 <- utils::packageVersion("SeuratObject") >= "5.0.0"
  expected_access_api <- if (seurat_object_v5) "layer" else "slot"
  expected_access_strategy <- if (seurat_object_v5) {
    "SeuratObject >= 5: layer first, slot fallback"
  } else {
    "SeuratObject < 5 or unavailable: slot first, layer fallback"
  }
  expect_identical(
    cache$effective_parameters$assay_access_api,
    expected_access_api
  )
  expect_identical(
    cache$effective_parameters$assay_access_strategy,
    expected_access_strategy
  )
})

test_that("CORALembedding rejects fewer than five training cells", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("irlba")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("matrixStats")

  object <- make_test_iscoral_object(n_cells = 20L, n_genes = 10L)
  cells <- colnames(object)
  probability <- Matrix::Diagonal(n = length(cells), x = 1)
  dimnames(probability) <- list(cells, cells)
  gene_statistics <- data.frame(
    omega_square_infer = seq(1, 0.1, length.out = nrow(object)),
    row.names = rownames(object)
  )
  object@misc$CORAL_scan_results <- list(
    binsize_200 = list(
      S_stab = probability,
      ss_filter_average = gene_statistics
    )
  )

  expect_error(
    CORALembedding(
      object,
      n_dims = 2,
      top_n_genes = 8,
      max_cells_lasso = 4,
      workers_glmnet = 1,
      random_seed = 47
    ),
    "at least 5 training cells"
  )
})

test_that("CORALembedding requires at least two matching memory genes", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  object <- make_test_iscoral_object(n_cells = 20L, n_genes = 10L)
  cells <- colnames(object)
  probability <- Matrix::Diagonal(n = length(cells), x = 1)
  dimnames(probability) <- list(cells, cells)
  gene_statistics <- data.frame(
    omega_square_infer = seq(1, 0.1, length.out = nrow(object)),
    row.names = rownames(object)
  )
  object@misc$CORAL_scan_results <- list(
    binsize_200 = list(
      S_stab = probability,
      ss_filter_average = gene_statistics
    )
  )

  expect_error(
    CORALembedding(object, n_dims = 2, top_n_genes = 1),
    "At least two ranked memory genes"
  )
})

test_that("CORALembedding validates both probability-matrix axes", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  object <- make_test_iscoral_object(n_cells = 20L, n_genes = 10L)
  cells <- colnames(object)
  probability <- Matrix::sparseMatrix(
    i = seq_along(cells),
    j = seq_along(cells),
    x = 1,
    dims = c(length(cells), length(cells)),
    dimnames = list(
      cells,
      c(cells[-length(cells)], "unknown_cell")
    )
  )
  gene_statistics <- data.frame(
    omega_square_infer = seq(1, 0.1, length.out = nrow(object)),
    row.names = rownames(object)
  )
  object@misc$CORAL_scan_results <- list(
    binsize_200 = list(
      S_stab = probability,
      ss_filter_average = gene_statistics
    )
  )

  expect_error(
    CORALembedding(object, n_dims = 2, top_n_genes = 8),
    "row and column cell-barcode sets"
  )
})
