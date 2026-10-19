#' Assign deterministic pseudoclones from the isCORAL consensus matrix
#'
#' @description
#' Converts the isCORAL consensus similarity matrix into deterministic
#' hierarchical pseudoclone assignments. The function aligns `S_stab` to the
#' current Seurat cell order, clusters the dissimilarity matrix `1 - S_stab`,
#' cuts the resulting tree at each requested value of `k`, and writes the
#' assignments to Seurat metadata.
#'
#' @param seurat_obj A Seurat object containing isCORAL scan results.
#' @param bin_size The bin size used in the isCORAL scan.
#' @param k_values A numeric vector of positive integer values used to cut the
#'   hierarchical tree.
#' @param hclust_method Agglomeration method passed to `stats::hclust()`.
#' @param ... Arguments forwarded by the deprecated `SimulateMetaClone()` alias.
#'
#' @return A modified Seurat object with metadata columns named
#'   `pseudoclone_k<k>_bs<bin_size>`.
#' @details
#' This function does not use random sampling. Given the same `S_stab`, cell
#' order, `k_values`, and `hclust_method`, it returns the same assignments. It
#' does not read or simulate gene-expression values; use
#' `simulate_clonal_expression()` to visualize observed expression across the
#' resulting pseudoclones. `SimulateClone()`,
#' `Simulate_clonal_expression_simple()`, and `SimulateMetaClone()` are retained
#' only as deprecated compatibility aliases.
#' @examples
#' \dontrun{
#' seurat_obj <- Assign_deterministic_pseudoclone(
#'   seurat_obj = seurat_obj,
#'   bin_size = 200,
#'   k_values = c(5, 10)
#' )
#' }
#' @export
Assign_deterministic_pseudoclone <- function(
    seurat_obj,
    bin_size,
    k_values,
    hclust_method = "ward.D2"
) {
  if (!"CORAL_scan_results" %in% names(seurat_obj@misc)) {
    stop("Error: isCORAL results not found in seurat_obj@misc$CORAL_scan_results.")
  }

  if (length(bin_size) != 1L || is.na(bin_size) || !nzchar(as.character(bin_size))) {
    stop("'bin_size' must contain one non-missing value.", call. = FALSE)
  }
  if (!is.numeric(k_values) || length(k_values) == 0L || anyNA(k_values) ||
      any(!is.finite(k_values)) || any(k_values < 1) ||
      any(k_values != floor(k_values))) {
    stop("'k_values' must contain positive integers.", call. = FALSE)
  }
  if (!is.character(hclust_method) || length(hclust_method) != 1L ||
      is.na(hclust_method) || !nzchar(hclust_method)) {
    stop("'hclust_method' must be one non-empty character value.", call. = FALSE)
  }

  result_name <- paste0("binsize_", bin_size)
  if (!result_name %in% names(seurat_obj@misc$CORAL_scan_results)) {
    stop(paste("Error: No isCORAL results found for bin_size =", bin_size))
  }

  coral_result <- seurat_obj@misc$CORAL_scan_results[[result_name]]
  if (!"S_stab" %in% names(coral_result)) {
    stop(paste("Error: 'S_stab' consensus matrix not found for bin_size =", bin_size))
  }

  message("Computing hierarchical clustering from consensus matrix (S_stab)...")
  S_mat <- as.matrix(coral_result$S_stab)
  seurat_cells <- colnames(seurat_obj)

  if (is.null(seurat_cells) || anyNA(seurat_cells) ||
      any(!nzchar(seurat_cells)) || anyDuplicated(seurat_cells)) {
    stop("Error: Seurat cell names must be non-empty and unique.")
  }
  if (length(seurat_cells) != nrow(S_mat) || length(seurat_cells) != ncol(S_mat)) {
    stop("Error: S_stab dimensions do not match the Seurat object.")
  }
  if (any(k_values > length(seurat_cells))) {
    stop("Every value in 'k_values' must be no larger than the number of cells.", call. = FALSE)
  }

  # Fill a missing axis only when dimensions match the current Seurat object.
  if (is.null(rownames(S_mat)) && is.null(colnames(S_mat))) {
    rownames(S_mat) <- seurat_cells
    colnames(S_mat) <- seurat_cells
  } else if (is.null(rownames(S_mat))) {
    rownames(S_mat) <- colnames(S_mat)
  } else if (is.null(colnames(S_mat))) {
    colnames(S_mat) <- rownames(S_mat)
  }

  if (anyNA(rownames(S_mat)) || anyNA(colnames(S_mat)) ||
      any(!nzchar(rownames(S_mat))) || any(!nzchar(colnames(S_mat))) ||
      anyDuplicated(rownames(S_mat)) || anyDuplicated(colnames(S_mat)) ||
      !setequal(rownames(S_mat), seurat_cells) ||
      !setequal(colnames(S_mat), seurat_cells)) {
    stop("Error: S_stab cell names do not match the Seurat object.")
  }

  # Reordering both axes makes the deterministic cell order explicit.
  S_mat <- S_mat[seurat_cells, seurat_cells, drop = FALSE]
  dist_matrix <- stats::as.dist(1 - S_mat)
  hclust_obj <- stats::hclust(dist_matrix, method = hclust_method)
  k_values <- unique(as.integer(k_values))

  message(paste(
    "Generating deterministic pseudoclone assignments from bin_size =",
    bin_size,
    "at k =",
    paste(k_values, collapse = ", ")
  ))

  for (k in k_values) {
    clusters <- as.factor(stats::cutree(hclust_obj, k = k))
    names(clusters) <- rownames(S_mat)

    col_name <- paste0("pseudoclone_k", k, "_bs", bin_size)
    seurat_obj[[col_name]] <- clusters
  }

  seurat_obj
}

#' Visualize gene expression across deterministic isCORAL pseudoclones
#'
#' @description
#' Reproduces the isCORAL-specific portion of the clone-expression benchmark
#' visualization for user-selected genes. For each gene, the function draws its
#' observed expression on an existing UMAP or other two-dimensional reduction
#' and a violin plot across deterministic isCORAL pseudoclones. It does not draw
#' experimental ground truth or comparator methods and does not generate
#' synthetic expression values.
#'
#' @param seurat_obj A Seurat object containing deterministic pseudoclone labels.
#' @param features Character vector of genes to visualize.
#' @param group_by Metadata column containing labels generated by
#'   `Assign_deterministic_pseudoclone()`.
#' @param assay Assay containing the expression values. If `NULL`, the current
#'   default assay is used.
#' @param layer Normalized-expression layer for Seurat v5, or slot for Seurat
#'   v4. The default is `"data"`.
#' @param reduction Two-dimensional reduction used for the expression map. If
#'   `NULL`, `"umap.raw"` is preferred over `"umap"`.
#' @param min_cells Minimum number of cells required for a pseudoclone to be
#'   retained.
#' @param n_pseudoclones Maximum number of pseudoclones displayed, selected by
#'   decreasing cell count. Use `NULL` to retain all eligible pseudoclones.
#' @param invalid_pseudoclones Labels excluded before pseudoclone-size filtering.
#' @param background_color Low-expression color for the embedding.
#' @param signal_color High-expression color and isCORAL focal color.
#' @param violin_colors Two colors defining the pseudoclone violin palette.
#' @param point_size Point size in the expression embedding.
#' @param base_size Base font size for the violin panel.
#' @param base_family Base font family.
#' @param return_data If `TRUE`, return the plot together with cell-level data,
#'   pseudoclone summaries, selected labels, and resolved parameters.
#'
#' @return A patchwork object, or a list when `return_data = TRUE`.
#' @details
#' Eligible pseudoclones are selected once by cell count so every requested gene
#' is shown on the same pseudoclone subset. Within each gene panel, pseudoclones
#' are ordered by decreasing mean expression, matching the benchmark workflow.
#' The expression map contains all cells shared by the Seurat object, selected
#' assay layer, metadata, and reduction; the violin panel contains only cells in
#' the selected pseudoclones.
#' @examples
#' \dontrun{
#' group_column <- "pseudoclone_k10_bs200"
#' p <- simulate_clonal_expression(
#'   seurat_obj = seurat_obj,
#'   features = c("NNMT", "HIST1H1E"),
#'   group_by = group_column,
#'   assay = "RNA",
#'   layer = "data",
#'   reduction = "umap"
#' )
#' print(p)
#' }
#' @importFrom rlang .data
#' @export
simulate_clonal_expression <- function(
    seurat_obj,
    features,
    group_by,
    assay = NULL,
    layer = "data",
    reduction = NULL,
    min_cells = 2,
    n_pseudoclones = 12,
    invalid_pseudoclones = c("NA", "NaN", "", "None", "Unknown", "unknown", "Unassigned"),
    background_color = "#E0E4E7",
    signal_color = "#C74B46",
    violin_colors = c("#F5C6C5", "#C74B46"),
    point_size = 0.8,
    base_size = 14,
    base_family = "Arial",
    return_data = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw clonal-expression plots.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required to compose clonal-expression plots.", call. = FALSE)
  }

  features <- unique(as.character(features))
  if (length(features) == 0L || anyNA(features) || any(!nzchar(features))) {
    stop("'features' must contain at least one non-empty gene name.", call. = FALSE)
  }
  if (!is.character(group_by) || length(group_by) != 1L ||
      is.na(group_by) || !nzchar(group_by)) {
    stop("'group_by' must be one non-empty metadata column name.", call. = FALSE)
  }
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Metadata column '%s' was not found.", group_by), call. = FALSE)
  }
  if (!is.character(layer) || length(layer) != 1L || is.na(layer) || !nzchar(layer)) {
    stop("'layer' must be one non-empty character value.", call. = FALSE)
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1L || is.na(min_cells) ||
      !is.finite(min_cells) || min_cells < 1 || min_cells != floor(min_cells)) {
    stop("'min_cells' must be one positive integer.", call. = FALSE)
  }
  if (!is.null(n_pseudoclones) &&
      (!is.numeric(n_pseudoclones) || length(n_pseudoclones) != 1L ||
       is.na(n_pseudoclones) || !is.finite(n_pseudoclones) ||
       n_pseudoclones < 1 || n_pseudoclones != floor(n_pseudoclones))) {
    stop("'n_pseudoclones' must be NULL or one positive integer.", call. = FALSE)
  }
  if (!is.character(violin_colors) || length(violin_colors) != 2L ||
      anyNA(violin_colors) || any(!nzchar(violin_colors))) {
    stop("'violin_colors' must contain exactly two non-empty colors.", call. = FALSE)
  }
  for (color_name in c("background_color", "signal_color")) {
    color_value <- get(color_name)
    if (!is.character(color_value) || length(color_value) != 1L ||
        is.na(color_value) || !nzchar(color_value)) {
      stop(sprintf("'%s' must be one non-empty color value.", color_name), call. = FALSE)
    }
  }
  if (!is.numeric(point_size) || length(point_size) != 1L ||
      is.na(point_size) || !is.finite(point_size) || point_size <= 0) {
    stop("'point_size' must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(base_size) || length(base_size) != 1L ||
      is.na(base_size) || !is.finite(base_size) || base_size <= 0) {
    stop("'base_size' must be one positive finite number.", call. = FALSE)
  }
  if (!is.character(base_family) || length(base_family) != 1L || is.na(base_family)) {
    stop("'base_family' must be one character value.", call. = FALSE)
  }
  if (!is.logical(return_data) || length(return_data) != 1L || is.na(return_data)) {
    stop("'return_data' must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seurat_obj)
  }
  if (!is.character(assay) || length(assay) != 1L || is.na(assay) || !nzchar(assay)) {
    stop("'assay' must be NULL or one non-empty assay name.", call. = FALSE)
  }
  if (!assay %in% names(seurat_obj@assays)) {
    stop(sprintf("Assay '%s' was not found.", assay), call. = FALSE)
  }

  expression_matrix <- .isCORAL_get_assay_layer(
    seurat_obj = seurat_obj,
    assay = assay,
    layer = layer,
    required = TRUE
  )
  missing_features <- setdiff(features, rownames(expression_matrix))
  if (length(missing_features) > 0L) {
    stop(sprintf(
      "Feature(s) not found in assay '%s' layer '%s': %s.",
      assay,
      layer,
      paste(missing_features, collapse = ", ")
    ), call. = FALSE)
  }

  reduction <- .iscoral_resolve_clonal_expression_reduction(seurat_obj, reduction)
  embedding <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  if (ncol(embedding) < 2L) {
    stop(sprintf("Reduction '%s' must contain at least two dimensions.", reduction), call. = FALSE)
  }

  plot_cells <- Reduce(
    intersect,
    list(
      colnames(seurat_obj),
      colnames(expression_matrix),
      rownames(seurat_obj@meta.data),
      rownames(embedding)
    )
  )
  if (length(plot_cells) < 2L) {
    stop("Fewer than two cells overlap across the Seurat object, expression layer, metadata, and reduction.", call. = FALSE)
  }

  expression_subset <- as.matrix(expression_matrix[features, plot_cells, drop = FALSE])
  if (any(!is.finite(expression_subset))) {
    stop("Requested expression values contain non-finite entries.", call. = FALSE)
  }
  embedding <- embedding[plot_cells, 1:2, drop = FALSE]
  group_labels <- as.character(seurat_obj@meta.data[plot_cells, group_by])

  valid_group <- !is.na(group_labels) & nzchar(group_labels) &
    !(group_labels %in% invalid_pseudoclones)
  pseudoclone_counts <- table(group_labels[valid_group])
  pseudoclone_counts <- pseudoclone_counts[pseudoclone_counts >= as.integer(min_cells)]
  if (length(pseudoclone_counts) == 0L) {
    stop("No pseudoclones passed the requested label and cell-count filters.", call. = FALSE)
  }

  count_table <- data.frame(
    Pseudoclone = names(pseudoclone_counts),
    Cells = as.integer(pseudoclone_counts),
    stringsAsFactors = FALSE
  )
  count_table <- count_table[
    order(-count_table$Cells, count_table$Pseudoclone),
    ,
    drop = FALSE
  ]
  if (!is.null(n_pseudoclones)) {
    count_table <- utils::head(count_table, as.integer(n_pseudoclones))
  }
  selected_pseudoclones <- count_table$Pseudoclone

  all_expression_data <- vector("list", length(features))
  violin_data <- vector("list", length(features))
  summary_data <- vector("list", length(features))
  pseudoclone_orders <- vector("list", length(features))
  plot_list <- vector("list", length(features))
  names(all_expression_data) <- features
  names(violin_data) <- features
  names(summary_data) <- features
  names(pseudoclone_orders) <- features
  names(plot_list) <- features

  for (feature_index in seq_along(features)) {
    feature_name <- features[[feature_index]]
    feature_values <- as.numeric(expression_subset[feature_name, ])

    feature_data <- data.frame(
      Cell = plot_cells,
      Feature = feature_name,
      Expression = feature_values,
      Pseudoclone = group_labels,
      Selected = valid_group & group_labels %in% selected_pseudoclones,
      stringsAsFactors = FALSE
    )
    all_expression_data[[feature_index]] <- feature_data

    feature_violin_data <- feature_data[feature_data$Selected, , drop = FALSE]
    feature_summary <- do.call(rbind, lapply(selected_pseudoclones, function(label) {
      values <- feature_violin_data$Expression[feature_violin_data$Pseudoclone == label]
      data.frame(
        Feature = feature_name,
        Pseudoclone = label,
        Cells = length(values),
        MeanExpression = mean(values),
        MedianExpression = stats::median(values),
        stringsAsFactors = FALSE
      )
    }))
    feature_summary <- feature_summary[
      order(-feature_summary$MeanExpression, feature_summary$Pseudoclone),
      ,
      drop = FALSE
    ]
    pseudoclone_order <- feature_summary$Pseudoclone
    feature_violin_data$Pseudoclone <- factor(
      feature_violin_data$Pseudoclone,
      levels = pseudoclone_order
    )

    violin_data[[feature_index]] <- feature_violin_data
    summary_data[[feature_index]] <- feature_summary
    pseudoclone_orders[[feature_index]] <- pseudoclone_order

    embedding_data <- data.frame(
      Dim1 = embedding[, 1],
      Dim2 = embedding[, 2],
      Expression = feature_values,
      Cell = plot_cells,
      stringsAsFactors = FALSE
    )
    # Draw low expression first so high-expression cells remain visible.
    embedding_data <- embedding_data[
      order(embedding_data$Expression, embedding_data$Cell),
      ,
      drop = FALSE
    ]

    p_feature <- ggplot2::ggplot(
      embedding_data,
      ggplot2::aes(x = .data$Dim1, y = .data$Dim2)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$Expression),
        size = point_size
      ) +
      ggplot2::scale_color_gradient(
        low = background_color,
        high = signal_color,
        name = "Expression"
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = paste(feature_name, "expression")) +
      ggplot2::theme_void(base_family = base_family) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = base_size + 2,
          hjust = 0.5,
          face = "plain"
        ),
        plot.subtitle = ggplot2::element_blank(),
        plot.caption = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(size = max(base_size - 3, 8)),
        legend.text = ggplot2::element_text(size = max(base_size - 4, 8)),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )

    pseudoclone_colors <- grDevices::colorRampPalette(violin_colors)(
      length(pseudoclone_order)
    )
    names(pseudoclone_colors) <- pseudoclone_order

    p_violin <- ggplot2::ggplot(
      feature_violin_data,
      ggplot2::aes(
        x = .data$Pseudoclone,
        y = .data$Expression,
        fill = .data$Pseudoclone
      )
    ) +
      ggplot2::geom_violin(
        scale = "width",
        trim = TRUE,
        alpha = 0.9,
        color = "#333333",
        linewidth = 0.5
      ) +
      ggplot2::scale_fill_manual(values = pseudoclone_colors, drop = FALSE) +
      ggplot2::labs(
        title = "isCORAL pseudoclones",
        x = "Pseudoclone",
        y = "Expression"
      ) +
      ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(
          size = base_size + 2,
          hjust = 0.5,
          face = "plain"
        ),
        plot.subtitle = ggplot2::element_blank(),
        plot.caption = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          color = "#333333",
          linewidth = 1.2,
          fill = NA
        ),
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          color = "black",
          size = max(base_size - 5, 8)
        ),
        axis.text.y = ggplot2::element_text(
          color = "black",
          size = max(base_size - 3, 8)
        ),
        axis.title = ggplot2::element_text(size = max(base_size - 1, 8)),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )

    plot_list[[feature_index]] <- patchwork::wrap_plots(
      p_feature,
      p_violin,
      nrow = 1,
      widths = c(1, 1.4)
    )
  }

  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 1)
  expression_data <- do.call(rbind, all_expression_data)
  violin_plot_data <- do.call(rbind, violin_data)
  pseudoclone_summary <- do.call(rbind, summary_data)
  rownames(expression_data) <- NULL
  rownames(violin_plot_data) <- NULL
  rownames(pseudoclone_summary) <- NULL

  if (return_data) {
    return(list(
      plot = combined_plot,
      plots = plot_list,
      expression_data = expression_data,
      violin_data = violin_plot_data,
      pseudoclone_summary = pseudoclone_summary,
      selected_pseudoclones = selected_pseudoclones,
      pseudoclone_order = pseudoclone_orders,
      group_by = group_by,
      assay = assay,
      layer = layer,
      reduction = reduction,
      min_cells = as.integer(min_cells),
      n_pseudoclones = if (is.null(n_pseudoclones)) NULL else as.integer(n_pseudoclones)
    ))
  }

  combined_plot
}

.iscoral_resolve_clonal_expression_reduction <- function(seurat_obj, reduction = NULL) {
  available_reductions <- names(seurat_obj@reductions)
  if (!is.null(reduction)) {
    if (!is.character(reduction) || length(reduction) != 1L ||
        is.na(reduction) || !nzchar(reduction)) {
      stop("'reduction' must be NULL or one non-empty reduction name.", call. = FALSE)
    }
    if (!reduction %in% available_reductions) {
      stop(sprintf(
        "Reduction '%s' was not found. Available reductions: %s.",
        reduction,
        paste(available_reductions, collapse = ", ")
      ), call. = FALSE)
    }
    return(reduction)
  }

  candidates <- c("umap.raw", "umap")
  available_candidates <- candidates[candidates %in% available_reductions]
  if (length(available_candidates) == 0L) {
    stop(
      "No original UMAP reduction was found. Please provide 'reduction' explicitly.",
      call. = FALSE
    )
  }
  available_candidates[[1]]
}
