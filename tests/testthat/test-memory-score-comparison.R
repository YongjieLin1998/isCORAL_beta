make_memory_score_comparison_object <- function() {
  counts <- matrix(
    1,
    nrow = 4,
    ncol = 4,
    dimnames = list(
      c("geneA", "geneB", "geneC", "geneD"),
      paste0("cell_", 1:4)
    )
  )
  counts["geneC", ] <- 0
  object <- Seurat::CreateSeuratObject(counts = counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object@misc$CORAL_scan_results <- list(
    binsize_200 = list(
      ss_filter_average = data.frame(
        name = c("geneA", "geneB", "geneC", "geneD"),
        mean = c(0.2, 0.3, 0.05, 0.4),
        omega_square_infer = c(0.2, 0.8, 0.9, NA_real_),
        stringsAsFactors = FALSE
      )
    )
  )
  object@misc$groundtruth <- list(
    heritable_genes_df = data.frame(
      name = c("geneA", "geneB", "geneC", "geneD"),
      Omega_square = c(0.1, 0.7, 0.8, 0.4),
      stringsAsFactors = FALSE
    ),
    heritable_genes_threshold = 0.5
  )
  object
}

test_that("CompareMemoryScores aligns and summarizes memory scores", {
  skip_if_not_installed("Seurat")

  object <- make_memory_score_comparison_object()
  result <- CompareMemoryScores(
    object,
    binsize_key = "binsize_200",
    assay = "RNA",
    min_mean = 0.1,
    top_n = 1
  )

  expect_true("CompareMemoryScores" %in% getNamespaceExports("isCORAL"))
  expect_identical(result$comparison$name, c("geneB", "geneA"))
  expect_identical(
    colnames(result$comparison),
    c(
      "name",
      "mean_expression",
      "inferred_memory_score",
      "ground_truth_memory_score",
      "ground_truth_rank",
      "is_high_truth",
      "above_permutation_threshold"
    )
  )
  expect_identical(result$comparison$ground_truth_rank, c(1L, 3L))
  expect_identical(result$comparison$is_high_truth, c(TRUE, FALSE))
  expect_identical(
    result$comparison$above_permutation_threshold,
    c(TRUE, FALSE)
  )
  expect_equal(result$summary$genes_compared, 2)
  expect_equal(result$summary$eligible_ground_truth_genes, 3)
  expect_equal(result$summary$positive_genes, 1)
  expect_equal(result$summary$negative_genes, 1)
  expect_equal(result$summary$min_mean, 0.1)
  expect_equal(result$summary$top_n, 1)
  expect_equal(result$summary$effective_top_n, 1)
  expect_equal(result$summary$spearman_r, 1)
  expect_equal(result$summary$pearson_r, 1)
  expect_equal(result$summary$rmse, 0.1)
  expect_equal(result$summary$auroc, 1)
  expect_equal(result$summary$auprc, 1)
  expect_equal(result$summary$random_auprc, 0.5)
  expect_equal(result$summary$auprc_enrichment, 2)
  expect_equal(result$summary$max_f1, 1)
  expect_equal(result$summary$best_threshold, 0.8)
  expect_equal(result$summary$ground_truth_threshold, 0.5)
  expect_equal(
    result$roc_curve$false_positive_rate,
    c(0, 0, 1)
  )
  expect_equal(result$roc_curve$true_positive_rate, c(0, 1, 1))
  expect_equal(result$pr_curve$recall, c(0, 1, 1))
  expect_equal(result$pr_curve$precision, c(1, 1, 0.5))
  expect_identical(result$binsize_key, "binsize_200")
  expect_identical(result$misc_key, "CORAL_scan_results")
  expect_identical(result$assay, "RNA")
  expect_identical(result$layer, "data")
  expect_identical(result$top_n, 1L)
  expect_s3_class(result$plots$scatter, "ggplot")
  expect_s3_class(result$plots$roc, "ggplot")
  expect_s3_class(result$plots$pr, "ggplot")
  expect_identical(result$plots$roc$data, result$roc_curve)
  expect_identical(result$plots$pr$data, result$pr_curve)
  expect_identical(result$plot_data$predicted_high, c(TRUE, FALSE))
  expect_identical(result$plot_data$classification, c("TP", "TN"))
  expect_identical(
    result$plot_data[, colnames(result$comparison), drop = FALSE],
    result$comparison
  )
  expect_equal(result$plot_parameters$benchmark_truth_cutoff, 0.7)
  expect_false(result$plot_parameters$truth_cutoff_tied)
})

test_that("CompareMemoryScores requires calculated experimental ground truth", {
  skip_if_not_installed("Seurat")

  object <- make_memory_score_comparison_object()
  object@misc$groundtruth <- NULL

  expect_error(
    CompareMemoryScores(object, binsize_key = "binsize_200"),
    "Run calculate_heritable_genes\\(\\) first"
  )
})

test_that("CompareMemoryScores validates the benchmark positive-set size", {
  skip_if_not_installed("Seurat")

  object <- make_memory_score_comparison_object()

  expect_error(
    CompareMemoryScores(
      object,
      binsize_key = "binsize_200",
      top_n = 0
    ),
    "'top_n' must be one positive integer"
  )
})

test_that("CompareMemoryScores reports an undefined single-class benchmark", {
  skip_if_not_installed("Seurat")

  object <- make_memory_score_comparison_object()
  expect_warning(
    result <- CompareMemoryScores(
      object,
      binsize_key = "binsize_200",
      top_n = 2000
    ),
    "contains only one benchmark class"
  )

  expect_true(is.na(result$summary$auroc))
  expect_true(is.na(result$summary$auprc))
  expect_true(is.na(result$summary$max_f1))
  expect_equal(result$summary$spearman_r, 1)
  expect_s3_class(result$plots$scatter, "ggplot")
  expect_null(result$plots$roc)
  expect_null(result$plots$pr)
  expect_true(all(is.na(result$plot_data$classification)))
  expect_true(all(is.na(result$plot_data$predicted_high)))
})

test_that("CompareMemoryScores follows the benchmark curve integration", {
  skip_if_not_installed("Seurat")

  object <- make_memory_score_comparison_object()
  score_table <- object@misc$CORAL_scan_results$binsize_200$ss_filter_average
  score_table$omega_square_infer <- c(0.8, 0.2, 0.9, NA_real_)
  object@misc$CORAL_scan_results$binsize_200$ss_filter_average <- score_table
  result <- CompareMemoryScores(
    object,
    binsize_key = "binsize_200",
    assay = "RNA",
    top_n = 1
  )

  expect_equal(result$summary$auroc, 0)
  expect_equal(result$summary$auprc, 0.75)
  expect_equal(result$summary$max_f1, 2 / 3)
  expect_equal(result$summary$best_threshold, 0.2)
})

test_that("benchmark plotting options preserve the numerical comparison", {
  skip_if_not_installed("Seurat")

  object <- make_memory_score_comparison_object()
  original_misc <- object@misc
  plotted <- CompareMemoryScores(object, top_n = 1, label_top_n = 0)
  tables_only <- CompareMemoryScores(object, top_n = 1, make_plots = FALSE)
  original_fields <- c(
    "comparison", "summary", "roc_curve", "pr_curve", "binsize_key",
    "misc_key", "assay", "layer", "min_mean", "top_n"
  )

  expect_identical(plotted[original_fields], tables_only[original_fields])
  expect_identical(plotted$plot_data, tables_only$plot_data)
  expect_identical(tables_only$plots, list(scatter = NULL, roc = NULL, pr = NULL))
  expect_identical(object@misc, original_misc)
  expect_error(CompareMemoryScores(object, make_plots = NA), "'make_plots'")
  expect_error(CompareMemoryScores(object, label_top_n = -1), "'label_top_n'")
  expect_error(CompareMemoryScores(object, base_size = 0), "'base_size'")
  expect_error(CompareMemoryScores(object, base_family = NA_character_), "'base_family'")
})

test_that("scatter classification preserves rank-based ground-truth cutoff ties", {
  skip_if_not_installed("Seurat")

  object <- make_memory_score_comparison_object()
  object@misc$groundtruth$heritable_genes_df$Omega_square <- c(0.7, 0.7, 0.8, 0.4)
  result <- CompareMemoryScores(object, top_n = 1)

  expect_identical(result$comparison$name, c("geneA", "geneB"))
  expect_identical(result$comparison$is_high_truth, c(TRUE, FALSE))
  expect_equal(result$plot_parameters$benchmark_truth_cutoff, 0.7)
  expect_true(result$plot_parameters$truth_cutoff_tied)
  expect_identical(result$plot_data$classification, c("TP", "FP"))
})
