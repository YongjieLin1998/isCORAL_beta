#' Visualize ranked isCORAL memory genes
#'
#' @description
#' Reads gene statistics from an isCORAL scan result, retains genes whose mean
#' expression is strictly greater than `min_mean`, ranks them by inferred
#' omega-squared memory score, and draws a horizontal bar chart. The stored
#' scan result is not modified.
#'
#' @param seurat_obj A Seurat object containing isCORAL scan results.
#' @param binsize_key Result key under `seurat_obj@misc$CORAL_scan_results`,
#'   for example `"binsize_200"`. If `NULL`, the function uses a top-level
#'   `ss_filter_average` table or the only nested result containing one.
#' @param min_mean Strict lower bound for the `mean` column. The default keeps
#'   genes with `mean > 0.1`.
#' @param top_n Number of highest-ranked genes to display. Use `NULL` to display
#'   all genes that pass `min_mean`.
#' @param bar_fill Bar fill color.
#' @param bar_color Bar outline color.
#' @param base_size Base font size passed to `ggplot2::theme_classic()`.
#' @param base_family Font family passed to `ggplot2::theme_classic()`.
#' @param title Plot title. Use `NULL` to omit the title.
#' @param return_data If `TRUE`, return a list containing the plot, complete
#'   filtered table, displayed table, and resolved scan-result information.
#'
#' @return A ggplot object, or a list when `return_data = TRUE`.
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' result <- PlotMemoryGenes(
#'   seurat_obj,
#'   binsize_key = "binsize_200",
#'   min_mean = 0.1,
#'   top_n = 30,
#'   return_data = TRUE
#' )
#' print(result$plot)
#' head(result$memory_genes)
#' }
#' @export
PlotMemoryGenes <- function(
    seurat_obj,
    binsize_key = NULL,
    min_mean = 0.1,
    top_n = 30,
    bar_fill = "#C74B46",
    bar_color = "#333333",
    base_size = 15,
    base_family = "Arial",
    title = "Top-ranked memory genes",
    return_data = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw memory-gene plots.", call. = FALSE)
  }
  if (!is.numeric(min_mean) || length(min_mean) != 1L ||
      is.na(min_mean) || !is.finite(min_mean)) {
    stop("'min_mean' must be one finite numeric value.", call. = FALSE)
  }
  if (!is.null(top_n) &&
      (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
       !is.finite(top_n) || top_n < 1 || top_n != floor(top_n))) {
    stop("'top_n' must be NULL or one positive integer.", call. = FALSE)
  }
  if (!is.numeric(base_size) || length(base_size) != 1L ||
      is.na(base_size) || !is.finite(base_size) || base_size <= 0) {
    stop("'base_size' must be one positive finite number.", call. = FALSE)
  }
  for (color_name in c("bar_fill", "bar_color")) {
    color_value <- get(color_name)
    if (!is.character(color_value) || length(color_value) != 1L ||
        is.na(color_value) || !nzchar(color_value)) {
      stop(sprintf("'%s' must be one non-empty color value.", color_name), call. = FALSE)
    }
  }
  if (!is.character(base_family) || length(base_family) != 1L ||
      is.na(base_family)) {
    stop("'base_family' must be one character value.", call. = FALSE)
  }
  if (!is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("'title' must be NULL or one character value.", call. = FALSE)
  }
  if (!is.logical(return_data) || length(return_data) != 1L || is.na(return_data)) {
    stop("'return_data' must be TRUE or FALSE.", call. = FALSE)
  }

  memory_info <- .iscoral_get_memory_gene_stats(seurat_obj, binsize_key)
  memory_gene_stats <- as.data.frame(memory_info$statistics)
  required_columns <- c("name", "mean", "omega_square_infer")
  missing_columns <- setdiff(required_columns, colnames(memory_gene_stats))
  if (length(missing_columns) > 0L) {
    stop(sprintf(
      "The memory-gene table is missing required column(s): %s.",
      paste(missing_columns, collapse = ", ")
    ), call. = FALSE)
  }
  if (!is.numeric(memory_gene_stats$mean) ||
      !is.numeric(memory_gene_stats$omega_square_infer)) {
    stop("The 'mean' and 'omega_square_infer' columns must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(memory_gene_stats$mean)) ||
      any(!is.finite(memory_gene_stats$omega_square_infer))) {
    stop("Memory-gene statistics contain non-finite values.", call. = FALSE)
  }

  memory_gene_stats$name <- as.character(memory_gene_stats$name)
  if (anyNA(memory_gene_stats$name) || any(!nzchar(memory_gene_stats$name)) ||
      anyDuplicated(memory_gene_stats$name)) {
    stop("The memory-gene table must contain unique, non-empty gene names.", call. = FALSE)
  }

  memory_genes <- memory_gene_stats[
    memory_gene_stats$mean > min_mean,
    ,
    drop = FALSE
  ]
  if (nrow(memory_genes) == 0L) {
    stop(sprintf("No memory genes passed mean > %s.", min_mean), call. = FALSE)
  }
  memory_genes <- memory_genes[
    order(
      -memory_genes$omega_square_infer,
      -memory_genes$mean,
      memory_genes$name
    ),
    ,
    drop = FALSE
  ]
  rownames(memory_genes) <- NULL

  n_display <- if (is.null(top_n)) nrow(memory_genes) else min(as.integer(top_n), nrow(memory_genes))
  plot_data <- utils::head(memory_genes, n_display)
  plot_data$gene <- factor(
    plot_data$name,
    levels = rev(plot_data$name)
  )

  memory_gene_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$gene, y = .data$omega_square_infer)
  ) +
    ggplot2::geom_col(
      width = 0.75,
      fill = bar_fill,
      color = bar_color,
      linewidth = 0.4
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = expression(omega^2 ~ memory ~ score)
    ) +
    ggplot2::theme_classic(
      base_size = base_size,
      base_family = base_family
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = base_size + 3,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        size = max(base_size - 1, 8)
      ),
      axis.text.y = ggplot2::element_text(
        size = max(base_size - 4, 8),
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(
        size = max(base_size - 3, 8),
        color = "black"
      ),
      axis.line = ggplot2::element_line(color = bar_color)
    )

  if (return_data) {
    return(list(
      plot = memory_gene_plot,
      memory_genes = memory_genes,
      plot_data = plot_data,
      binsize_key = memory_info$binsize_key,
      misc_key = memory_info$misc_key,
      min_mean = min_mean,
      top_n = top_n
    ))
  }

  memory_gene_plot
}

#' Compare inferred memory scores with experimental lineage ground truth
#'
#' @description
#' Aligns isCORAL-inferred memory scores with the gene-wise omega-squared values
#' previously calculated from experimental lineage labels by
#' `calculate_heritable_genes()`. The comparison is restricted to genes shared
#' by both results and the selected log-normalized `data` layer whose mean
#' expression is strictly greater than `min_mean`.
#'
#' @param seurat_obj A Seurat object containing both isCORAL scan results and
#'   `seurat_obj@misc$groundtruth$heritable_genes_df`.
#' @param binsize_key Result key under `seurat_obj@misc$CORAL_scan_results`,
#'   for example `"binsize_200"`. If `NULL`, the same unambiguous result
#'   resolution used by `PlotMemoryGenes()` is applied.
#' @param assay Assay used to calculate gene means from its log-normalized
#'   `data` layer. If `NULL`, the current default assay is used.
#' @param min_mean Strict lower bound applied to the mean expression calculated
#'   from the selected assay's `data` layer. The default reproduces the
#'   benchmark filter `mean > 0.1`.
#' @param top_n Number of highest-ranking experimental omega-squared genes used
#'   as the positive class for the benchmark metrics. The default reproduces
#'   the benchmark definition of the top 2,000 memory genes.
#'
#' @return A list containing `comparison`, the aligned gene-level table;
#'   `summary`, the benchmark metrics and resolved parameters; `roc_curve` and
#'   `pr_curve`, the threshold-level curve coordinates; and the resolved scan
#'   result identifiers.
#' @details
#' This function does not calculate experimental ground truth and does not
#' modify the Seurat object. Run `calculate_heritable_genes()` first with the
#' experimental lineage-label metadata column. Spearman correlation, Pearson
#' correlation, and RMSE compare the continuous inferred and ground-truth
#' omega-squared scores. AUROC, AUPRC, AUPRC enrichment, and maximum F1 use the
#' top experimental omega-squared genes as the positive class, following the
#' package's memory-gene benchmark workflow. If the filtered comparison has no
#' positive or no negative genes, classification metrics are returned as `NA`
#' with a warning, while continuous-score metrics remain available.
#' @examples
#' \dontrun{
#' seurat_obj <- calculate_heritable_genes(
#'   seurat_obj,
#'   group_by = "barcode",
#'   permutation_repeats = 100
#' )
#' comparison <- CompareMemoryScores(
#'   seurat_obj,
#'   binsize_key = "binsize_200",
#'   assay = "RNA",
#'   top_n = 2000
#' )
#' comparison$summary
#' head(comparison$comparison)
#' }
#' @export
CompareMemoryScores <- function(
    seurat_obj,
    binsize_key = NULL,
    assay = NULL,
    min_mean = 0.1,
    top_n = 2000
) {
  if (!is.numeric(min_mean) || length(min_mean) != 1L ||
      is.na(min_mean) || !is.finite(min_mean)) {
    stop("'min_mean' must be one finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
      !is.finite(top_n) || top_n < 1 || top_n != floor(top_n) ||
      top_n > .Machine$integer.max) {
    stop("'top_n' must be one positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seurat_obj)
  }
  if (!is.character(assay) || length(assay) != 1L ||
      is.na(assay) || !nzchar(assay)) {
    stop("'assay' must be NULL or one non-empty assay name.", call. = FALSE)
  }
  if (!assay %in% names(seurat_obj@assays)) {
    stop(sprintf("Assay '%s' was not found.", assay), call. = FALSE)
  }

  memory_info <- .iscoral_get_memory_gene_stats(seurat_obj, binsize_key)
  inferred_table <- as.data.frame(memory_info$statistics)
  inferred_columns <- c("name", "omega_square_infer")
  missing_inferred <- setdiff(inferred_columns, colnames(inferred_table))
  if (length(missing_inferred) > 0L) {
    stop(sprintf(
      "The inferred memory-gene table is missing required column(s): %s.",
      paste(missing_inferred, collapse = ", ")
    ), call. = FALSE)
  }

  groundtruth <- seurat_obj@misc$groundtruth
  if (is.null(groundtruth) || is.null(groundtruth$heritable_genes_df)) {
    stop(
      paste0(
        "Experimental memory scores were not found in ",
        "seurat_obj@misc$groundtruth$heritable_genes_df. ",
        "Run calculate_heritable_genes() first."
      ),
      call. = FALSE
    )
  }
  groundtruth_table <- as.data.frame(groundtruth$heritable_genes_df)
  required_groundtruth_columns <- c("name", "Omega_square")
  missing_groundtruth <- setdiff(
    required_groundtruth_columns,
    colnames(groundtruth_table)
  )
  if (length(missing_groundtruth) > 0L) {
    stop(sprintf(
      "The ground-truth memory table is missing required column(s): %s.",
      paste(missing_groundtruth, collapse = ", ")
    ), call. = FALSE)
  }

  inferred_table$name <- as.character(inferred_table$name)
  groundtruth_table$name <- as.character(groundtruth_table$name)
  if (anyNA(inferred_table$name) || any(!nzchar(inferred_table$name)) ||
      anyDuplicated(inferred_table$name)) {
    stop(
      "The inferred memory-gene table must contain unique, non-empty gene names.",
      call. = FALSE
    )
  }
  if (anyNA(groundtruth_table$name) || any(!nzchar(groundtruth_table$name)) ||
      anyDuplicated(groundtruth_table$name)) {
    stop(
      "The ground-truth memory table must contain unique, non-empty gene names.",
      call. = FALSE
    )
  }
  if (!is.numeric(inferred_table$omega_square_infer) ||
      !is.numeric(groundtruth_table$Omega_square)) {
    stop(
      "The inferred memory score and ground-truth omega-squared columns must be numeric.",
      call. = FALSE
    )
  }
  groundtruth_table$.ground_truth_order <- seq_len(nrow(groundtruth_table))
  groundtruth_columns <- c(
    required_groundtruth_columns,
    ".ground_truth_order"
  )

  expression_matrix <- .isCORAL_get_assay_layer(
    seurat_obj = seurat_obj,
    assay = assay,
    layer = "data",
    required = TRUE
  )
  expression_means <- data.frame(
    name = rownames(expression_matrix),
    mean_expression = as.numeric(Matrix::rowMeans(expression_matrix)),
    stringsAsFactors = FALSE
  )
  if (anyNA(expression_means$name) || any(!nzchar(expression_means$name)) ||
      anyDuplicated(expression_means$name)) {
    stop(
      "The selected assay must contain unique, non-empty gene names.",
      call. = FALSE
    )
  }

  eligible_genes <- merge(
    expression_means,
    groundtruth_table[, groundtruth_columns, drop = FALSE],
    by = "name",
    all = FALSE,
    sort = FALSE
  )
  eligible_genes <- eligible_genes[
    is.finite(eligible_genes$mean_expression) &
      eligible_genes$mean_expression > min_mean &
      is.finite(eligible_genes$Omega_square),
    ,
    drop = FALSE
  ]
  if (nrow(eligible_genes) < 2L) {
    stop(
      "Fewer than two shared genes passed the requested mean and finite-value filters.",
      call. = FALSE
    )
  }

  eligible_genes <- eligible_genes[
    order(
      -eligible_genes$Omega_square,
      eligible_genes$.ground_truth_order
    ),
    ,
    drop = FALSE
  ]
  eligible_genes$ground_truth_rank <- seq_len(nrow(eligible_genes))
  effective_top_n <- min(top_n, nrow(eligible_genes))
  eligible_genes$is_high_truth <-
    eligible_genes$ground_truth_rank <= effective_top_n

  comparison <- merge(
    eligible_genes,
    inferred_table[, c("name", "omega_square_infer"), drop = FALSE],
    by = "name",
    all = FALSE,
    sort = FALSE
  )
  comparison <- comparison[
    is.finite(comparison$omega_square_infer),
    ,
    drop = FALSE
  ]
  if (nrow(comparison) < 2L) {
    stop(
      "Fewer than two shared genes have finite inferred memory scores.",
      call. = FALSE
    )
  }

  colnames(comparison)[colnames(comparison) == "omega_square_infer"] <-
    "inferred_memory_score"
  colnames(comparison)[colnames(comparison) == "Omega_square"] <-
    "ground_truth_memory_score"

  threshold <- groundtruth$heritable_genes_threshold
  if (is.null(threshold) || !is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || !is.finite(threshold)) {
    threshold <- NA_real_
  } else {
    threshold <- as.numeric(threshold)
  }
  comparison$above_permutation_threshold <- if (is.finite(threshold)) {
    comparison$ground_truth_memory_score > threshold
  } else {
    rep(NA, nrow(comparison))
  }
  comparison <- comparison[, c(
    "name",
    "mean_expression",
    "inferred_memory_score",
    "ground_truth_memory_score",
    "ground_truth_rank",
    "is_high_truth",
    "above_permutation_threshold"
  ), drop = FALSE]

  comparison <- comparison[
    order(
      -comparison$ground_truth_memory_score,
      comparison$ground_truth_rank,
      comparison$name
    ),
    ,
    drop = FALSE
  ]
  rownames(comparison) <- NULL

  safe_correlation <- function(x, y, method) {
    if (length(unique(x)) < 2L || length(unique(y)) < 2L) {
      return(NA_real_)
    }
    unname(stats::cor(x, y, method = method))
  }

  if (all(comparison$is_high_truth) || !any(comparison$is_high_truth)) {
    warning(
      paste0(
        "The filtered comparison contains only one benchmark class; ",
        "AUROC, AUPRC, AUPRC enrichment, maximum F1, and the best ",
        "threshold are returned as NA. Choose a smaller 'top_n' or review ",
        "the input filters to obtain both classes."
      ),
      call. = FALSE
    )
  }

  curve_results <- .iscoral_memory_score_curves(
    labels = comparison$is_high_truth,
    scores = comparison$inferred_memory_score
  )

  summary_table <- data.frame(
    genes_compared = nrow(comparison),
    eligible_ground_truth_genes = nrow(eligible_genes),
    positive_genes = sum(comparison$is_high_truth),
    negative_genes = sum(!comparison$is_high_truth),
    min_mean = min_mean,
    top_n = top_n,
    effective_top_n = effective_top_n,
    assay = assay,
    layer = "data",
    spearman_r = safe_correlation(
      comparison$inferred_memory_score,
      comparison$ground_truth_memory_score,
      method = "spearman"
    ),
    pearson_r = safe_correlation(
      comparison$inferred_memory_score,
      comparison$ground_truth_memory_score,
      method = "pearson"
    ),
    rmse = sqrt(mean(
      (comparison$inferred_memory_score -
         comparison$ground_truth_memory_score)^2
    )),
    auroc = curve_results$auroc,
    auprc = curve_results$auprc,
    random_auprc = curve_results$random_auprc,
    auprc_enrichment = curve_results$auprc_enrichment,
    max_f1 = curve_results$max_f1,
    best_threshold = curve_results$best_threshold,
    ground_truth_threshold = threshold,
    stringsAsFactors = FALSE
  )

  list(
    comparison = comparison,
    summary = summary_table,
    roc_curve = curve_results$roc_curve,
    pr_curve = curve_results$pr_curve,
    binsize_key = memory_info$binsize_key,
    misc_key = memory_info$misc_key,
    assay = assay,
    layer = "data",
    min_mean = min_mean,
    top_n = top_n
  )
}

.iscoral_memory_score_curves <- function(labels, scores) {
  labels <- as.logical(labels)
  scores <- as.numeric(scores)
  n_positive <- sum(labels)
  n_negative <- sum(!labels)
  random_auprc <- n_positive / length(labels)

  if (n_positive == 0L || n_negative == 0L) {
    return(list(
      roc_curve = data.frame(
        false_positive_rate = c(0, 1),
        true_positive_rate = c(0, 1),
        threshold = c(Inf, -Inf)
      ),
      pr_curve = data.frame(
        recall = c(0, 1),
        precision = c(1, random_auprc),
        threshold = c(Inf, -Inf),
        f1 = c(0, 2 * random_auprc / (1 + random_auprc))
      ),
      auroc = NA_real_,
      auprc = NA_real_,
      random_auprc = random_auprc,
      auprc_enrichment = NA_real_,
      max_f1 = NA_real_,
      best_threshold = NA_real_
    ))
  }

  score_order <- order(scores, decreasing = TRUE)
  sorted_scores <- scores[score_order]
  sorted_labels <- labels[score_order]
  group_ends <- which(!duplicated(sorted_scores, fromLast = TRUE))
  true_positives <- c(0, cumsum(sorted_labels)[group_ends])
  false_positives <- c(0, cumsum(!sorted_labels)[group_ends])
  thresholds <- c(Inf, sorted_scores[group_ends])

  true_positive_rate <- true_positives / n_positive
  false_positive_rate <- false_positives / n_negative
  auroc <- sum(
    diff(false_positive_rate) *
      (true_positive_rate[-1] +
         true_positive_rate[-length(true_positive_rate)]) / 2
  )
  roc_curve <- data.frame(
    false_positive_rate = false_positive_rate,
    true_positive_rate = true_positive_rate,
    threshold = thresholds
  )

  recall <- true_positive_rate
  precision <- c(
    1,
    true_positives[-1] /
      (true_positives[-1] + false_positives[-1])
  )
  unique_recall <- !duplicated(recall)
  recall_for_auc <- recall[unique_recall]
  precision_for_auc <- precision[unique_recall]
  auprc <- sum(
    diff(recall_for_auc) *
      (precision_for_auc[-1] +
         precision_for_auc[-length(precision_for_auc)]) / 2
  )
  f1 <- ifelse(
    precision + recall == 0,
    0,
    2 * precision * recall / (precision + recall)
  )
  best_index <- which.max(f1)
  pr_curve <- data.frame(
    recall = recall,
    precision = precision,
    threshold = thresholds,
    f1 = f1
  )

  list(
    roc_curve = roc_curve,
    pr_curve = pr_curve,
    auroc = auroc,
    auprc = auprc,
    random_auprc = random_auprc,
    auprc_enrichment = auprc / random_auprc,
    max_f1 = f1[[best_index]],
    best_threshold = thresholds[[best_index]]
  )
}

.iscoral_get_memory_gene_stats <- function(seurat_obj, binsize_key = NULL) {
  misc_names <- names(seurat_obj@misc)
  misc_key <- NULL
  if ("CORAL_scan_results" %in% misc_names) {
    misc_key <- "CORAL_scan_results"
  } else if ("coral_scan_results" %in% misc_names) {
    misc_key <- "coral_scan_results"
  }
  if (is.null(misc_key)) {
    stop(
      "Could not find isCORAL scan results in seurat_obj@misc$CORAL_scan_results.",
      call. = FALSE
    )
  }

  scan_results <- seurat_obj@misc[[misc_key]]
  if (!is.null(binsize_key)) {
    if (!is.character(binsize_key) || length(binsize_key) != 1L ||
        is.na(binsize_key) || !nzchar(binsize_key)) {
      stop("'binsize_key' must be NULL or one non-empty character value.", call. = FALSE)
    }
    statistics <- scan_results[[binsize_key]]$ss_filter_average
    if (is.null(statistics)) {
      stop(sprintf(
        "Could not find ss_filter_average under seurat_obj@misc$%s[['%s']].",
        misc_key,
        binsize_key
      ), call. = FALSE)
    }
    return(list(
      statistics = statistics,
      binsize_key = binsize_key,
      misc_key = misc_key
    ))
  }

  if (!is.null(scan_results$ss_filter_average)) {
    return(list(
      statistics = scan_results$ss_filter_average,
      binsize_key = NULL,
      misc_key = misc_key
    ))
  }

  result_keys <- names(scan_results)[vapply(scan_results, function(result) {
    is.list(result) && !is.null(result$ss_filter_average)
  }, logical(1))]
  if (length(result_keys) == 1L) {
    resolved_key <- result_keys[[1]]
    return(list(
      statistics = scan_results[[resolved_key]]$ss_filter_average,
      binsize_key = resolved_key,
      misc_key = misc_key
    ))
  }
  if (length(result_keys) > 1L) {
    stop(sprintf(
      paste0(
        "Multiple isCORAL result keys contain ss_filter_average (%s). ",
        "Please specify 'binsize_key'."
      ),
      paste(result_keys, collapse = ", ")
    ), call. = FALSE)
  }

  stop(sprintf(
    "Could not find any ss_filter_average table under seurat_obj@misc$%s.",
    misc_key
  ), call. = FALSE)
}
