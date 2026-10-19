test_that("PlotMemoryGenes filters and ranks memory genes", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  counts <- matrix(
    c(3, 4, 5, 6, 1, 1, 2, 1, 5, 5, 4, 6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("gene_a", "gene_b", "gene_c"),
      paste0("cell_", 1:4)
    )
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  object@misc$CORAL_scan_results <- list(
    binsize_200 = list(
      ss_filter_average = data.frame(
        mean = c(0.2, 0.05, 0.3),
        sst = c(1, 1, 1),
        name = c("gene_a", "gene_b", "gene_c"),
        fano = c(1, 1, 1),
        omega_square_infer = c(0.3, 0.9, 0.7),
        stringsAsFactors = FALSE
      )
    )
  )

  result <- PlotMemoryGenes(
    object,
    binsize_key = "binsize_200",
    min_mean = 0.1,
    top_n = 1,
    return_data = TRUE
  )

  expect_s3_class(result$plot, "ggplot")
  expect_identical(result$memory_genes$name, c("gene_c", "gene_a"))
  expect_identical(result$plot_data$name, "gene_c")
  expect_identical(result$binsize_key, "binsize_200")
  expect_identical(result$misc_key, "CORAL_scan_results")
  expect_identical(result$min_mean, 0.1)
  expect_identical(result$top_n, 1)
})

test_that("PlotMemoryGenes requires an explicit key for multiple scan results", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  counts <- matrix(
    1,
    nrow = 1,
    ncol = 3,
    dimnames = list("gene_a", paste0("cell_", 1:3))
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  statistics <- data.frame(
    mean = 0.2,
    name = "gene_a",
    omega_square_infer = 0.5,
    stringsAsFactors = FALSE
  )
  object@misc$CORAL_scan_results <- list(
    binsize_100 = list(ss_filter_average = statistics),
    binsize_200 = list(ss_filter_average = statistics)
  )

  expect_error(
    PlotMemoryGenes(object),
    "Multiple isCORAL result keys"
  )
})
