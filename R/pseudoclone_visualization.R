#' Visualize isCORAL pseudoclone probability fields
#'
#' @description
#' Draws a compact UMAP or CORAL-embedding visualization of one or two
#' isCORAL pseudoclones from the consensus matrix `S_stab`. With one seed cell,
#' the function shows that seed's probability/enrichment field. With two seed
#' cells, it shows a red/blue winner-takes-most competition map. If `clone_col`
#' is supplied, an optional ground-truth clone panel is placed beside the
#' isCORAL prediction.
#'
#' @param seurat_obj A Seurat object containing isCORAL scan results.
#' @param seed_cell Optional seed cell barcode, or a character vector of one or
#'   two seed cell barcodes.
#' @param seed_cells Deprecated plural alias for `seed_cell`. If both are
#'   `NULL`, one seed cell is sampled at random.
#' @param clone_col Optional metadata column containing true clone labels. Use
#'   `NULL` when ground truth is unavailable.
#' @param clone_ids Optional character vector of one or two clone IDs used to
#'   choose seed cells from `clone_col`.
#' @param binsize_key Result key under `seurat_obj@misc$CORAL_scan_results`,
#'   for example `"binsize_200"`. No default is assumed; if this is `NULL`,
#'   the function uses a top-level `S_stab` or the only available nested result.
#' @param prob_mat Optional precomputed cell-cell probability/consensus matrix.
#'   When supplied, this matrix is used instead of reading `S_stab` from
#'   `seurat_obj@misc`.
#' @param reduction Dimensional reduction to plot. If `NULL`, the function uses
#'   the first available reduction among `"umap.raw"`, `"umap"`,
#'   `"CORAL_clone_embed"`, and `"CORAL_cell_embed"`.
#'   In clone space, background and inferred cells use open circles;
#'   seed and ground-truth markers retain their filled shapes.
#' @param invalid_clones Clone labels to ignore when selecting seeds by clone.
#' @param min_enrichment Cells below this fold-enrichment are shown as
#'   background.
#' @param max_enrichment Upper limit for alpha scaling.
#' @param colors Two colors used for seed A and seed B.
#' @param background_color Color for background cells.
#' @param point_size Point size for regular cells.
#' @param seed_size Point size for seed cells.
#' @param show_ground_truth Whether to draw the clone ground-truth panel when
#'   `clone_col` is available.
#' @param title Title for the isCORAL prediction panel.
#' @param random_seed Optional random seed used only when seed cells are sampled.
#' @param seed Deprecated alias for `random_seed`.
#' @param return_data If `TRUE`, return a list containing the plot and plotting
#'   data instead of only the plot object.
#'
#' @return A ggplot/patchwork object, or a list when `return_data = TRUE`.
#' @importFrom rlang .data
#' @importFrom utils head
#' @examples
#' \dontrun{
#' p <- PlotPseudoclone(
#'   seurat_obj,
#'   clone_col = "true_clone",
#'   clone_ids = c("11", "27"),
#'   binsize_key = "binsize_200",
#'   reduction = "umap.raw"
#' )
#'
#' p <- PlotPseudoclone(
#'   seurat_obj,
#'   seed_cell = c("cell_A", "cell_B"),
#'   binsize_key = "binsize_200",
#'   reduction = "CORAL_clone_embed",
#'   show_ground_truth = FALSE
#' )
#'
#' p <- PlotPseudoclone(
#'   seurat_obj,
#'   reduction = "umap.raw",
#'   random_seed = 1
#' )
#' }
#' @export
PlotPseudoclone <- function(
    seurat_obj,
    seed_cell = NULL,
    seed_cells = NULL,
    clone_col = NULL,
    clone_ids = NULL,
    binsize_key = NULL,
    prob_mat = NULL,
    reduction = NULL,
    invalid_clones = c("0", "NA", "NaN", "", "None", "Unknown", "unknown", "Unassigned"),
    min_enrichment = 5,
    max_enrichment = 200,
    colors = c("#C74B46", "#3D5A80"),
    background_color = "#E0E0E0",
    point_size = 1.2,
    seed_size = 5,
    show_ground_truth = !is.null(clone_col),
    title = "isCORAL pseudoclone",
    random_seed = NULL,
    seed = NULL,
    return_data = FALSE
) {
  clone_col <- .iscoral_normalize_clone_col(clone_col)
  show_ground_truth <- isTRUE(show_ground_truth) && !is.null(clone_col)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw pseudoclone plots.", call. = FALSE)
  }
  if (show_ground_truth && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required when show_ground_truth = TRUE.", call. = FALSE)
  }
  if (length(colors) < 2) {
    stop("'colors' must contain at least two colors.", call. = FALSE)
  }
  if (!is.numeric(min_enrichment) || length(min_enrichment) != 1 || min_enrichment <= 0) {
    stop("'min_enrichment' must be one positive number.", call. = FALSE)
  }
  if (!is.numeric(max_enrichment) || length(max_enrichment) != 1 || max_enrichment <= min_enrichment) {
    stop("'max_enrichment' must be larger than 'min_enrichment'.", call. = FALSE)
  }

  if (!is.null(clone_col) && !clone_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("clone_col '%s' was not found in seurat_obj@meta.data.", clone_col), call. = FALSE)
  }
  if (!is.null(seed_cell) && !is.null(seed_cells)) {
    stop("Use only one of 'seed_cell' and 'seed_cells'.", call. = FALSE)
  }
  if (!is.null(seed_cells)) {
    warning("'seed_cells' is deprecated; use 'seed_cell'.", call. = FALSE)
    seed_cell <- seed_cells
  }
  if (!is.null(seed_cell)) {
    seed_cells <- seed_cell
  }
  if (!is.null(seed) && !is.null(random_seed)) {
    stop("Use only one of 'seed' and 'random_seed'.", call. = FALSE)
  }
  if (!is.null(seed)) {
    warning("'seed' is deprecated; use 'random_seed'.", call. = FALSE)
    random_seed <- seed
  }

  reduction <- .iscoral_choose_reduction(seurat_obj, reduction)
  point_shape <- if (identical(reduction, "CORAL_clone_embed")) 1 else 19
  embedding <- Seurat::Embeddings(seurat_obj, reduction)
  if (ncol(embedding) < 2) {
    stop(sprintf("Reduction '%s' must contain at least two dimensions.", reduction), call. = FALSE)
  }

  matrix_info <- .iscoral_get_pseudoclone_matrix(seurat_obj, binsize_key, prob_mat)
  prob_mat <- matrix_info$matrix
  binsize_key <- matrix_info$binsize_key
  prob_mat <- .iscoral_align_pseudoclone_matrix(prob_mat, seurat_obj)
  prob_mat <- .iscoral_normalize_probability_rows(prob_mat)

  plot_cells <- Reduce(
    intersect,
    list(rownames(embedding), rownames(prob_mat), colnames(prob_mat), colnames(seurat_obj))
  )
  if (length(plot_cells) < 2) {
    stop("No overlapping cells were found among the reduction, Seurat object, and probability matrix.", call. = FALSE)
  }

  embedding <- embedding[plot_cells, 1:2, drop = FALSE]
  prob_mat <- prob_mat[plot_cells, plot_cells, drop = FALSE]

  seed_cells <- .iscoral_resolve_seed_cells(
    seurat_obj = seurat_obj,
    seed_cells = seed_cells,
    clone_col = clone_col,
    clone_ids = clone_ids,
    valid_cells = plot_cells,
    invalid_clones = invalid_clones,
    seed = random_seed
  )

  n_seeds <- length(seed_cells)
  if (!n_seeds %in% c(1, 2)) {
    stop("Please provide one or two seed cells/clones.", call. = FALSE)
  }

  seed_clones <- rep(NA_character_, n_seeds)
  if (!is.null(clone_col)) {
    seed_clones <- as.character(seurat_obj@meta.data[seed_cells, clone_col])
  }

  plot_df <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Cell_ID = rownames(embedding),
    stringsAsFactors = FALSE
  )

  seed_df <- data.frame(
    Dim1 = embedding[seed_cells, 1],
    Dim2 = embedding[seed_cells, 2],
    Seed = factor(paste0("Seed ", seq_len(n_seeds)), levels = paste0("Seed ", seq_len(n_seeds))),
    Cell_ID = seed_cells,
    stringsAsFactors = FALSE
  )

  enrichments <- lapply(seed_cells, function(cell) {
    values <- as.numeric(prob_mat[cell, ]) * (ncol(prob_mat) - 1)
    names(values) <- colnames(prob_mat)
    values[cell] <- NA_real_
    values
  })

  prediction_df <- .iscoral_build_prediction_df(
    plot_df = plot_df,
    enrichments = enrichments,
    min_enrichment = min_enrichment
  )

  p_pred <- .iscoral_plot_prediction(
    prediction_df = prediction_df,
    seed_df = seed_df,
    n_seeds = n_seeds,
    colors = colors,
    background_color = background_color,
    point_size = point_size,
    seed_size = seed_size,
    min_enrichment = min_enrichment,
    max_enrichment = max_enrichment,
    title = title,
    point_shape = point_shape
  )

  ground_truth_df <- NULL
  plot_out <- p_pred
  if (show_ground_truth && !is.null(clone_col)) {
    ground_truth_df <- .iscoral_build_ground_truth_df(
      plot_df = plot_df,
      seurat_obj = seurat_obj,
      clone_col = clone_col,
      seed_clones = seed_clones,
      invalid_clones = invalid_clones
    )
    p_gt <- .iscoral_plot_ground_truth(
      ground_truth_df = ground_truth_df,
      seed_df = seed_df,
      seed_clones = seed_clones,
      n_seeds = n_seeds,
      colors = colors,
      background_color = background_color,
      point_size = point_size,
      seed_size = seed_size,
      point_shape = point_shape
    )
    plot_out <- patchwork::wrap_plots(p_gt, p_pred, nrow = 1, guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
  }

  if (return_data) {
    return(list(
      plot = plot_out,
      prediction_data = prediction_df,
      ground_truth_data = ground_truth_df,
      seed_cells = seed_cells,
      seed_clones = seed_clones,
      reduction = reduction,
      binsize_key = binsize_key
    ))
  }

  plot_out
}

.iscoral_choose_reduction <- function(seurat_obj, reduction = NULL) {
  available <- names(seurat_obj@reductions)
  if (!is.null(reduction)) {
    if (!reduction %in% available) {
      stop(sprintf("Reduction '%s' was not found. Available reductions: %s",
                   reduction, paste(available, collapse = ", ")), call. = FALSE)
    }
    return(reduction)
  }

  candidates <- c("umap.raw", "umap", "CORAL_clone_embed", "CORAL_cell_embed")
  chosen <- candidates[candidates %in% available]
  if (length(chosen) == 0) {
    stop("No default plotting reduction found. Please provide 'reduction'.", call. = FALSE)
  }
  chosen[[1]]
}

.iscoral_normalize_clone_col <- function(clone_col) {
  if (is.null(clone_col) || length(clone_col) == 0 || is.na(clone_col[[1]])) {
    return(NULL)
  }
  clone_col <- as.character(clone_col[[1]])
  if (clone_col %in% c("", "NULL", "NA", "None", "none")) {
    return(NULL)
  }
  clone_col
}

.iscoral_get_pseudoclone_matrix <- function(seurat_obj, binsize_key, prob_mat = NULL) {
  if (!is.null(prob_mat)) {
    return(list(matrix = prob_mat, binsize_key = binsize_key))
  }

  misc_names <- names(seurat_obj@misc)
  misc_key <- NULL
  if ("CORAL_scan_results" %in% misc_names) {
    misc_key <- "CORAL_scan_results"
  } else if ("coral_scan_results" %in% misc_names) {
    misc_key <- "coral_scan_results"
  }

  if (is.null(misc_key)) {
    stop("Could not find isCORAL scan results in seurat_obj@misc$CORAL_scan_results.", call. = FALSE)
  }

  scan_results <- seurat_obj@misc[[misc_key]]
  if (!is.null(binsize_key)) {
    if (!is.null(scan_results[[binsize_key]]$S_stab)) {
      return(list(matrix = scan_results[[binsize_key]]$S_stab, binsize_key = binsize_key))
    }
    stop(sprintf("Could not find S_stab under seurat_obj@misc$%s[['%s']].", misc_key, binsize_key), call. = FALSE)
  }

  if (!is.null(scan_results$S_stab)) {
    return(list(matrix = scan_results$S_stab, binsize_key = NULL))
  }

  result_keys <- names(scan_results)[vapply(scan_results, function(x) {
    is.list(x) && !is.null(x$S_stab)
  }, logical(1))]

  if (length(result_keys) == 1) {
    return(list(matrix = scan_results[[result_keys[[1]]]]$S_stab, binsize_key = result_keys[[1]]))
  }
  if (length(result_keys) > 1) {
    stop(sprintf(
      "Multiple isCORAL result keys contain S_stab (%s). Please specify 'binsize_key'.",
      paste(result_keys, collapse = ", ")
    ), call. = FALSE)
  }

  stop(sprintf("Could not find any S_stab under seurat_obj@misc$%s.", misc_key), call. = FALSE)
}

.iscoral_align_pseudoclone_matrix <- function(prob_mat, seurat_obj) {
  cell_names <- colnames(seurat_obj)
  matrix_has_names <- !is.null(rownames(prob_mat)) && !is.null(colnames(prob_mat))
  numeric_rows <- matrix_has_names && all(grepl("^[0-9]+$", head(rownames(prob_mat), 50)))

  if (!matrix_has_names || numeric_rows || length(intersect(rownames(prob_mat), cell_names)) == 0) {
    if (nrow(prob_mat) == length(cell_names) && ncol(prob_mat) == length(cell_names)) {
      rownames(prob_mat) <- cell_names
      colnames(prob_mat) <- cell_names
    } else {
      stop("The probability matrix has no usable cell names and its dimensions do not match the Seurat object.", call. = FALSE)
    }
  }

  common_cells <- Reduce(intersect, list(cell_names, rownames(prob_mat), colnames(prob_mat)))
  if (length(common_cells) < 2) {
    stop("Fewer than two cells overlap between the probability matrix and the Seurat object.", call. = FALSE)
  }

  prob_mat[common_cells, common_cells, drop = FALSE]
}

.iscoral_normalize_probability_rows <- function(prob_mat) {
  if (!inherits(prob_mat, "CsparseMatrix")) {
    prob_mat <- methods::as(prob_mat, "CsparseMatrix")
  }

  row_sums <- Matrix::rowSums(prob_mat)
  row_sums[row_sums == 0] <- 1
  prob_norm <- Matrix::Diagonal(x = 1 / row_sums) %*% prob_mat
  rownames(prob_norm) <- rownames(prob_mat)
  colnames(prob_norm) <- colnames(prob_mat)
  prob_norm
}

.iscoral_resolve_seed_cells <- function(
    seurat_obj,
    seed_cells = NULL,
    clone_col = NULL,
    clone_ids = NULL,
    valid_cells,
    invalid_clones,
    seed = NULL
) {
  if (!is.null(seed_cells)) {
    seed_cells <- as.character(seed_cells)
    missing_cells <- setdiff(seed_cells, valid_cells)
    if (length(missing_cells) > 0) {
      stop(sprintf("Seed cell(s) not present in the aligned plotting space: %s",
                   paste(missing_cells, collapse = ", ")), call. = FALSE)
    }
    return(seed_cells)
  }

  if (is.null(clone_col)) {
    if (!is.null(clone_ids)) {
      stop("Cannot use 'clone_ids' when 'clone_col' is NULL.", call. = FALSE)
    }
    return(.iscoral_sample_cells(valid_cells, n = 1, seed = seed))
  }

  meta <- seurat_obj@meta.data
  meta_cells <- intersect(rownames(meta), valid_cells)
  labels <- as.character(meta[meta_cells, clone_col])
  valid_label <- !is.na(labels) & !labels %in% invalid_clones

  if (!is.null(clone_ids)) {
    clone_ids <- as.character(clone_ids)
    if (!length(clone_ids) %in% c(1, 2)) {
      stop("'clone_ids' must contain one or two clone labels.", call. = FALSE)
    }
  } else {
    candidate_cells <- meta_cells[valid_label]
    if (length(candidate_cells) == 0) {
      candidate_cells <- valid_cells
    }
    return(.iscoral_sample_cells(candidate_cells, n = 1, seed = seed))
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  vapply(clone_ids, function(id) {
    candidates <- meta_cells[valid_label & labels == id]
    if (length(candidates) == 0) {
      stop(sprintf("Clone '%s' has no valid cells in the aligned plotting space.", id), call. = FALSE)
    }
    sample(candidates, 1)
  }, character(1))
}

.iscoral_sample_cells <- function(cells, n = 1, seed = NULL) {
  cells <- unique(as.character(cells))
  if (length(cells) < n) {
    stop(sprintf("Need at least %d valid cell(s), but only found %d.", n, length(cells)), call. = FALSE)
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  sample(cells, n)
}

.iscoral_build_prediction_df <- function(plot_df, enrichments, min_enrichment) {
  n_seeds <- length(enrichments)
  out <- plot_df

  if (n_seeds == 1) {
    out$Enrichment <- enrichments[[1]][match(out$Cell_ID, names(enrichments[[1]]))]
    out$Enrichment[is.na(out$Enrichment)] <- 0
    out$Assigned_Lineage <- ifelse(out$Enrichment >= min_enrichment, "Seed A", "Background")
    out$Dominant_E <- out$Enrichment
  } else {
    out$E_A <- enrichments[[1]][match(out$Cell_ID, names(enrichments[[1]]))]
    out$E_B <- enrichments[[2]][match(out$Cell_ID, names(enrichments[[2]]))]
    out$E_A[is.na(out$E_A)] <- 0
    out$E_B[is.na(out$E_B)] <- 0
    out$Dominant_E <- pmax(out$E_A, out$E_B)
    out$Assigned_Lineage <- ifelse(out$E_A >= out$E_B, "Seed A", "Seed B")
    out$Assigned_Lineage[out$Dominant_E < min_enrichment] <- "Background"
  }

  out[order(out$Dominant_E), , drop = FALSE]
}

.iscoral_plot_prediction <- function(
    prediction_df,
    seed_df,
    n_seeds,
    colors,
    background_color,
    point_size,
    seed_size,
    min_enrichment,
    max_enrichment,
    title,
    point_shape = 19
) {
  lineage_colors <- c("Seed A" = colors[[1]], "Seed B" = colors[[2]])
  seed_shapes <- c(21, 24)
  alpha_breaks <- unique(c(min_enrichment, 20, 100, max_enrichment))
  alpha_breaks <- alpha_breaks[alpha_breaks >= min_enrichment & alpha_breaks <= max_enrichment]

  p <- ggplot2::ggplot(prediction_df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2)) +
    ggplot2::geom_point(
      data = prediction_df[prediction_df$Assigned_Lineage == "Background", , drop = FALSE],
      color = background_color,
      size = point_size,
      shape = point_shape
    ) +
    ggplot2::geom_point(
      data = prediction_df[prediction_df$Assigned_Lineage != "Background", , drop = FALSE],
      ggplot2::aes(color = .data$Assigned_Lineage, alpha = .data$Dominant_E),
      size = point_size * 1.25,
      shape = point_shape
    ) +
    ggplot2::scale_color_manual(
      name = "Pseudoclone",
      values = lineage_colors[seq_len(n_seeds)]
    ) +
    ggplot2::scale_alpha_continuous(
      name = "Enrichment",
      limits = c(min_enrichment, max_enrichment),
      oob = .iscoral_oob_squish,
      range = c(0.15, 0.95),
      trans = "log1p",
      breaks = alpha_breaks
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1, size = 3.5)),
      alpha = ggplot2::guide_legend(order = 2, override.aes = list(color = "#555555", size = 3))
    ) +
    ggplot2::labs(title = title)

  for (i in seq_len(n_seeds)) {
    p <- p +
      ggplot2::geom_point(
        data = seed_df[i, , drop = FALSE],
        ggplot2::aes(x = .data$Dim1, y = .data$Dim2),
        shape = seed_shapes[[i]],
        size = seed_size,
        fill = colors[[i]],
        color = "black",
        stroke = 1.4,
        inherit.aes = FALSE
      )
  }

  p
}

.iscoral_build_ground_truth_df <- function(
    plot_df,
    seurat_obj,
    clone_col,
    seed_clones,
    invalid_clones
) {
  out <- plot_df
  labels <- as.character(seurat_obj@meta.data[out$Cell_ID, clone_col])
  labels[is.na(labels) | labels %in% invalid_clones] <- NA_character_

  out$GT_Label <- "Other"
  if (length(seed_clones) == 1) {
    out$GT_Label[!is.na(labels) & labels == seed_clones[[1]]] <- "Seed A clone"
  } else {
    out$GT_Label[!is.na(labels) & labels == seed_clones[[1]]] <- "Seed A clone"
    out$GT_Label[!is.na(labels) & labels == seed_clones[[2]]] <- "Seed B clone"
  }
  out
}

.iscoral_plot_ground_truth <- function(
    ground_truth_df,
    seed_df,
    seed_clones,
    n_seeds,
    colors,
    background_color,
    point_size,
    seed_size,
    point_shape = 19
) {
  gt_colors <- c("Seed A clone" = colors[[1]], "Seed B clone" = colors[[2]])
  seed_shapes <- c(21, 24)

  clone_text <- seed_clones[!is.na(seed_clones)]
  clone_text <- if (length(clone_text) == 0) "seed clones" else paste(clone_text, collapse = " & ")

  p <- ggplot2::ggplot(ground_truth_df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2)) +
    ggplot2::geom_point(
      data = ground_truth_df[ground_truth_df$GT_Label == "Other", , drop = FALSE],
      color = background_color,
      size = point_size,
      shape = point_shape
    ) +
    ggplot2::geom_point(
      data = ground_truth_df[ground_truth_df$GT_Label != "Other", , drop = FALSE],
      ggplot2::aes(fill = .data$GT_Label),
      shape = 21,
      size = point_size * 2.2,
      color = "black",
      stroke = 0.35
    ) +
    ggplot2::scale_fill_manual(name = "Ground truth", values = gt_colors[seq_len(n_seeds)]) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(title = sprintf("Ground truth\n%s", clone_text))

  for (i in seq_len(n_seeds)) {
    p <- p +
      ggplot2::geom_point(
        data = seed_df[i, , drop = FALSE],
        ggplot2::aes(x = .data$Dim1, y = .data$Dim2),
        shape = seed_shapes[[i]],
        size = seed_size,
        fill = colors[[i]],
        color = "black",
        stroke = 1.4,
        inherit.aes = FALSE
      )
  }

  p
}

.iscoral_oob_squish <- function(x, range = c(0, 1), ...) {
  pmin(pmax(x, range[[1]]), range[[2]])
}
