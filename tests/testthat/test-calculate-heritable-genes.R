make_heritable_gene_test_object <- function() {
  counts <- matrix(
    c(
      12, 11, 10, 1, 0, 1,
      10, 1, 9, 0, 8, 1,
      1, 2, 3, 8, 9, 10,
      2, 7, 2, 7, 2, 7
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("gene_", 1:4),
      paste0("cell_", 1:6)
    )
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object$barcode_a <- c("A", "A", "A", "B", "B", "B")
  object$barcode_b <- c("X", "X", "Y", "Y", "Z", "Z")
  object
}

test_that("calculate_heritable_genes never reuses stale ground-truth input", {
  skip_if_not_installed("Seurat")

  object <- make_heritable_gene_test_object()
  object@misc$groundtruth <- list(
    internal_data = list(
      sc_expression_subset = matrix(
        0,
        nrow = 1,
        ncol = 2,
        dimnames = list("stale_gene", c("stale_1", "stale_2"))
      ),
      barcode_used = c("stale", "stale")
    )
  )

  first_result <- calculate_heritable_genes(
    object,
    group_by = "barcode_a",
    permutation_repeats = 1,
    n_cores = 1,
    chunk_size = 2,
    verbose = FALSE
  )
  repeated_result <- calculate_heritable_genes(
    first_result,
    group_by = "barcode_b",
    invalid_barcodes = "Z",
    permutation_repeats = 1,
    n_cores = 1,
    chunk_size = 2,
    verbose = FALSE
  )
  fresh_result <- calculate_heritable_genes(
    make_heritable_gene_test_object(),
    group_by = "barcode_b",
    invalid_barcodes = "Z",
    permutation_repeats = 1,
    n_cores = 1,
    chunk_size = 2,
    verbose = FALSE
  )

  expect_null(first_result@misc$groundtruth$internal_data)
  expect_null(repeated_result@misc$groundtruth$internal_data)
  expect_identical(
    repeated_result@misc$groundtruth$input_provenance$group_by,
    "barcode_b"
  )
  expect_identical(
    repeated_result@misc$groundtruth$input_provenance$n_cells_input,
    6L
  )
  expect_identical(
    repeated_result@misc$groundtruth$input_provenance$n_cells_used,
    4L
  )
  expect_identical(
    repeated_result@misc$groundtruth$input_provenance$n_cells_removed,
    2L
  )
  expect_equal(
    repeated_result@misc$groundtruth$heritable_genes_df,
    fresh_result@misc$groundtruth$heritable_genes_df
  )
})
