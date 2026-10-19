make_test_hierarchical_clone_object <- function(n_cells = 12L, n_genes = 8L) {
  set.seed(211)
  counts <- matrix(
    stats::rpois(n_cells * n_genes, lambda = 4),
    nrow = n_genes,
    dimnames = list(
      paste0("gene_", seq_len(n_genes)),
      paste0("cell_", seq_len(n_cells))
    )
  )
  object <- Seurat::CreateSeuratObject(counts = counts)

  coordinates <- matrix(stats::rnorm(n_cells * 3L), nrow = n_cells)
  similarity <- exp(-as.matrix(stats::dist(coordinates)))
  dimnames(similarity) <- list(colnames(object), colnames(object))

  row_order <- rev(seq_len(n_cells))
  col_order <- c(seq.int(2L, n_cells), 1L)
  object@misc$CORAL_scan_results <- list(
    binsize_200 = list(
      S_stab = similarity[row_order, col_order, drop = FALSE]
    )
  )
  object
}

test_that("deterministic pseudoclone entry points remain compatible", {
  skip_if_not_installed("Seurat")

  expect_true(all(c(
    "Assign_deterministic_pseudoclone",
    "simulate_clonal_expression",
    "Simulate_clonal_expression_simple",
    "SimulateClone",
    "SimulateMetaClone"
  ) %in% getNamespaceExports("isCORAL")))
  expect_identical(
    formals(Assign_deterministic_pseudoclone),
    formals(SimulateClone)
  )

  object <- make_test_hierarchical_clone_object()
  canonical <- Assign_deterministic_pseudoclone(
    object,
    bin_size = 200,
    k_values = c(2, 3)
  )
  former_canonical <- expect_warning(
    Simulate_clonal_expression_simple(object, 200, c(2, 3)),
    "deprecated"
  )
  legacy_clone <- expect_warning(
    SimulateClone(object, 200, c(2, 3)),
    "deprecated"
  )
  legacy_meta <- expect_warning(
    SimulateMetaClone(object, 200, c(2, 3)),
    "deprecated"
  )

  expected_columns <- c("pseudoclone_k2_bs200", "pseudoclone_k3_bs200")
  expect_true(all(expected_columns %in% colnames(canonical[[]])))
  expect_false(any(grepl("^metaclone_k", colnames(canonical[[]]))))
  for (column_name in expected_columns) {
    expect_identical(
      as.character(canonical[[column_name]][[1]]),
      as.character(legacy_clone[[column_name]][[1]])
    )
    expect_identical(
      as.character(canonical[[column_name]][[1]]),
      as.character(former_canonical[[column_name]][[1]])
    )
    expect_identical(
      as.character(canonical[[column_name]][[1]]),
      as.character(legacy_meta[[column_name]][[1]])
    )
  }
  expect_identical(colnames(canonical), colnames(object))
})

test_that("deterministic pseudoclone assignment validates its cache", {
  skip_if_not_installed("Seurat")
  object <- make_test_hierarchical_clone_object()

  missing_scan <- object
  missing_scan@misc$CORAL_scan_results <- NULL
  expect_error(
    Assign_deterministic_pseudoclone(missing_scan, 200, 2),
    "CORAL_scan_results"
  )

  expect_error(
    Assign_deterministic_pseudoclone(object, 100, 2),
    "No isCORAL results"
  )

  missing_matrix <- object
  missing_matrix@misc$CORAL_scan_results$binsize_200$S_stab <- NULL
  expect_error(
    Assign_deterministic_pseudoclone(missing_matrix, 200, 2),
    "S_stab"
  )
})

test_that("clonal-expression API returns isCORAL-only plots and data", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("patchwork")

  object <- make_test_hierarchical_clone_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  set.seed(212)
  umap <- matrix(
    stats::rnorm(ncol(object) * 2L),
    ncol = 2L,
    dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))
  )
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = umap,
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )
  object <- Assign_deterministic_pseudoclone(object, 200, 3)

  result <- simulate_clonal_expression(
    seurat_obj = object,
    features = c("gene_1", "gene_2"),
    group_by = "pseudoclone_k3_bs200",
    assay = "RNA",
    layer = "data",
    reduction = "umap",
    min_cells = 2,
    n_pseudoclones = 3,
    return_data = TRUE
  )

  expect_s3_class(result$plot, "patchwork")
  expect_identical(names(result$plots), c("gene_1", "gene_2"))
  expect_true(all(c(
    "Cell", "Feature", "Expression", "Pseudoclone", "Selected"
  ) %in% colnames(result$expression_data)))
  expect_true(all(c(
    "Feature", "Pseudoclone", "Cells", "MeanExpression", "MedianExpression"
  ) %in% colnames(result$pseudoclone_summary)))
  expect_identical(result$group_by, "pseudoclone_k3_bs200")
  expect_identical(result$layer, "data")
})

test_that("clonal-expression API validates requested inputs", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("patchwork")

  object <- make_test_hierarchical_clone_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  set.seed(213)
  umap <- matrix(
    stats::rnorm(ncol(object) * 2L),
    ncol = 2L,
    dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))
  )
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = umap,
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(object)
  )
  object <- Assign_deterministic_pseudoclone(object, 200, 3)

  expect_error(
    simulate_clonal_expression(
      object,
      features = "missing_gene",
      group_by = "pseudoclone_k3_bs200"
    ),
    "not found"
  )
  expect_error(
    simulate_clonal_expression(
      object,
      features = "gene_1",
      group_by = "missing_column"
    ),
    "Metadata column"
  )
})
