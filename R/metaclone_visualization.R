#' Plot the seeds and inferred sisters of one Metaclone
#'
#' @description
#' Draws one Metaclone using the probability-cloud style of the embedding
#' scripts: all cells form a gray background, inferred sisters are colored
#' by probability-dependent opacity, and seed cells are drawn in a darker
#' shade of the same color.
#'
#' @param seurat_obj A Seurat object containing Metaclone metadata and the
#'   probability cache produced by `CORALembedding()`.
#' @param metaclone_id One Metaclone label from `cluster_col`.
#' @param reduction Reduction to plot. Defaults to `"CORAL_cell_embed"`.
#'   In `"CORAL_clone_embed"`, background and sister cells use open circles;
#'   seeds remain filled.
#' @param cluster_col Metadata column containing Metaclone labels.
#' @param dims Two embedding dimensions to display.
#' @param color Color of the inferred sister-cell cloud.
#' @param seed_darken Multiplier applied to the cloud color's RGB channels
#'   to obtain the darker seed color. Use 0.5 for the single-Metaclone
#'   embedding example or 0.45 for the two-Metaclone example.
#' @param background_color Color of the complete background cell layer.
#' @param background_size,point_size,seed_size Point sizes for background,
#'   inferred sisters, and seeds, respectively.
#' @param cell_type_col Optional metadata column used for density contours.
#'   Supply together with `contour_colors`.
#' @param contour_colors Named colors for the cell types to contour.
#' @param contour_level Normalized density contour level.
#' @param base_size,base_family Base font size and font family.
#' @param title Plot title. `NULL` uses the selected Metaclone label.
#' @param return_data If `TRUE`, return the plot, cell-level probabilities,
#'   seed and sister cell IDs, and plotting parameters.
#'
#' @details
#' All cells assigned to the selected Metaclone are seeds. The function
#' averages their rows of `seurat_obj@misc$isCORAL_embedding$P_norm`, sets
#' the probabilities at seed cells to zero, and divides by the maximum
#' remaining value. All non-seed cells with positive values contribute to
#' the inferred sister cloud. This is the calculation used in the original
#' `embedding/5.6.R` and `embedding/5.12.R` scripts.
#'
#' Probability rows, columns, metadata, and coordinates are matched by cell
#' name. An unnamed cache axis uses `cache$valid_cells`. Averaging and maximum
#' scaling use all cached target cells, including when the object has been
#' subsetted. Cells with coordinates but no cached probability remain gray
#' and retain `NA` probabilities in the returned table.
#'
#' The function uses the existing embedding and does not modify the Seurat
#' object, print a figure, or save a file.
#'
#' @return A ggplot object. With `return_data = TRUE`, a list containing
#'   `plot`, `plot_data`, `seed_cells`, `sister_cells`, `metaclone_id`,
#'   `reduction`, `cluster_col`, `source_binsize`, and `plot_parameters`.
#'   `plot_data` includes `MeanProbability` before seed exclusion,
#'   `SisterProbability` after seed exclusion, and `RelativeProbability`
#'   after maximum scaling.
#' @examples
#' \dontrun{
#' p <- PlotMetaclone(seurat_obj, metaclone_id = "2")
#' print(p)
#'
#' p <- PlotMetaclone(
#'   seurat_obj, metaclone_id = "2",
#'   cell_type_col = "cluster_integrated",
#'   contour_colors = c(HSC = "#D88DBA", Mk = "#F0AE46", GMP = "#3FA989")
#' )
#' print(p)
#' }
#' @export
PlotMetaclone <- function(
    seurat_obj,
    metaclone_id,
    reduction = "CORAL_cell_embed",
    cluster_col = "metaclone",
    dims = c(1, 2),
    color = "#3182BD",
    seed_darken = 0.5,
    background_color = "#E0E4E7",
    background_size = 0.8,
    point_size = 1.0,
    seed_size = 2.0,
    cell_type_col = NULL,
    contour_colors = NULL,
    contour_level = 0.15,
    base_size = 18,
    base_family = "Arial",
    title = NULL,
    return_data = FALSE
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.", call. = FALSE)
  }
  if (length(metaclone_id) != 1L || is.na(metaclone_id) ||
      !nzchar(as.character(metaclone_id))) {
    stop("Provide one non-missing 'metaclone_id'.", call. = FALSE)
  }
  metaclone_id <- as.character(metaclone_id)
  meta <- seurat_obj[[]]
  if (!is.character(cluster_col) || length(cluster_col) != 1L ||
      is.na(cluster_col) || !cluster_col %in% colnames(meta)) {
    stop("'cluster_col' must name an existing metadata column.", call. = FALSE)
  }
  if (!is.character(reduction) || length(reduction) != 1L ||
      is.na(reduction) || !reduction %in% names(seurat_obj@reductions)) {
    stop("The requested reduction was not found in the Seurat object.", call. = FALSE)
  }
  embedding <- Seurat::Embeddings(seurat_obj, reduction)
  if (!is.numeric(dims) || length(dims) != 2L || anyNA(dims) ||
      any(!is.finite(dims)) || any(dims != floor(dims)) ||
      any(dims < 1 | dims > ncol(embedding)) || anyDuplicated(dims)) {
    stop("'dims' must select two different embedding dimensions.", call. = FALSE)
  }
  positive_parameters <- list(
    background_size = background_size, point_size = point_size,
    seed_size = seed_size, base_size = base_size
  )
  for (parameter in names(positive_parameters)) {
    value <- positive_parameters[[parameter]]
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop(sprintf("'%s' must be one positive number.", parameter), call. = FALSE)
    }
  }
  if (!is.numeric(seed_darken) || length(seed_darken) != 1L ||
      !is.finite(seed_darken) || seed_darken < 0 || seed_darken > 1) {
    stop("'seed_darken' must be between 0 and 1.", call. = FALSE)
  }
  if (!is.character(color) || length(color) != 1L || is.na(color) ||
      !is.character(background_color) || length(background_color) != 1L ||
      is.na(background_color)) {
    stop("Provide one 'color' and one 'background_color'.", call. = FALSE)
  }
  if (!is.character(base_family) || length(base_family) != 1L || is.na(base_family)) {
    stop("'base_family' must be one font-family name.", call. = FALSE)
  }
  if (!is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("'title' must be NULL or one character string.", call. = FALSE)
  }
  if (!is.logical(return_data) || length(return_data) != 1L || is.na(return_data)) {
    stop("'return_data' must be TRUE or FALSE.", call. = FALSE)
  }

  result <- .iscoral_metaclone_plot_data(
    seurat_obj, metaclone_id, cluster_col, embedding, dims
  )
  plot_df <- result$plot_data

  # Match the darkened seed color used by the embedding examples.
  rgb_value <- as.numeric(grDevices::col2rgb(color)) * seed_darken
  seed_color <- grDevices::rgb(rgb_value[1], rgb_value[2], rgb_value[3], maxColorValue = 255)
  grDevices::col2rgb(background_color)

  cloud_df <- plot_df[plot_df$IsSister, , drop = FALSE]
  cloud_df <- cloud_df[order(cloud_df$RelativeProbability), , drop = FALSE]
  seed_df <- plot_df[plot_df$IsSeed, , drop = FALSE]
  point_shape <- if (identical(reduction, "CORAL_clone_embed")) 1 else 16
  if (is.null(title)) title <- paste("Metaclone", metaclone_id)

  # A full gray underlay preserves the manifold beneath the transparent cloud.
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = plot_df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2),
      color = background_color, size = background_size, shape = point_shape
    ) +
    ggplot2::geom_point(
      data = cloud_df,
      ggplot2::aes(x = .data$Dim1, y = .data$Dim2, alpha = .data$RelativeProbability),
      color = color, size = point_size, shape = point_shape
    ) +
    ggplot2::geom_point(
      data = seed_df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2),
      color = seed_color, size = seed_size, shape = 16
    ) +
    ggplot2::scale_alpha_continuous(range = c(0, 1), guide = "none") +
    ggplot2::theme_void(base_family = base_family, base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "plain", size = base_size,
        margin = ggplot2::margin(b = 10)
      ),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::labs(title = title)

  if (!is.null(cell_type_col) || !is.null(contour_colors)) {
    if (!is.character(cell_type_col) || length(cell_type_col) != 1L ||
        is.na(cell_type_col) || !cell_type_col %in% colnames(meta)) {
      stop("Provide an existing 'cell_type_col' together with 'contour_colors'.", call. = FALSE)
    }
    if (!is.character(contour_colors) || length(contour_colors) == 0L ||
        anyNA(contour_colors) || is.null(names(contour_colors)) ||
        anyNA(names(contour_colors)) || any(!nzchar(names(contour_colors))) ||
        anyDuplicated(names(contour_colors))) {
      stop("'contour_colors' must be a color vector named by cell type.", call. = FALSE)
    }
    if (!is.numeric(contour_level) || length(contour_level) != 1L ||
        !is.finite(contour_level) || contour_level <= 0 || contour_level >= 1) {
      stop("'contour_level' must be between 0 and 1, excluding the endpoints.", call. = FALSE)
    }
    grDevices::col2rgb(contour_colors)
    plot_df$CellType <- as.character(meta[plot_df$Cell_ID, cell_type_col])
    contour_df <- plot_df[plot_df$CellType %in% names(contour_colors), , drop = FALSE]
    for (cell_type in names(contour_colors)) {
      group <- contour_df[contour_df$CellType == cell_type, , drop = FALSE]
      if (nrow(group) < 3L || length(unique(group$Dim1)) < 2L ||
          length(unique(group$Dim2)) < 2L) {
        stop(sprintf("Cell type '%s' needs at least three cells and variation on both axes for a density contour.",
                     cell_type), call. = FALSE)
      }
    }
    p <- p +
      ggplot2::stat_density_2d(
        data = contour_df,
        ggplot2::aes(x = .data$Dim1, y = .data$Dim2, color = .data$CellType),
        geom = "path", contour_var = "ndensity", breaks = contour_level,
        linetype = "dashed", linewidth = 1.0, alpha = 0.9
      ) +
      ggplot2::scale_color_manual(values = contour_colors, guide = "none")
  }

  if (return_data) {
    return(list(
      plot = p,
      plot_data = plot_df,
      seed_cells = result$seed_cells,
      sister_cells = plot_df$Cell_ID[plot_df$IsSister],
      metaclone_id = metaclone_id,
      reduction = reduction,
      cluster_col = cluster_col,
      source_binsize = seurat_obj@misc$isCORAL_embedding$source_binsize,
      plot_parameters = list(
        dims = dims, color = color, seed_color = seed_color,
        seed_darken = seed_darken, background_color = background_color,
        background_size = background_size, point_size = point_size,
        seed_size = seed_size, cell_type_col = cell_type_col,
        contour_colors = contour_colors, contour_level = contour_level,
        base_size = base_size, base_family = base_family, title = title,
        probability_max = result$probability_max
      )
    ))
  }
  p
}

.iscoral_metaclone_plot_data <- function(
    seurat_obj, metaclone_id, cluster_col, embedding, dims
) {
  cache <- seurat_obj@misc$isCORAL_embedding
  prob_mat <- cache$P_norm
  if (is.null(prob_mat)) {
    stop("P_norm was not found; run CORALembedding() first.", call. = FALSE)
  }
  if (length(dim(prob_mat)) != 2L || nrow(prob_mat) != ncol(prob_mat)) {
    stop("Cached P_norm must be a square cell-cell matrix.", call. = FALSE)
  }
  valid_names <- function(x, n) {
    !is.null(x) && length(x) == n && !anyNA(x) &&
      all(nzchar(x)) && !anyDuplicated(x)
  }
  if (is.null(rownames(prob_mat)) || is.null(colnames(prob_mat))) {
    if (!valid_names(cache$valid_cells, nrow(prob_mat))) {
      stop("Unnamed P_norm axes require matching cache$valid_cells.", call. = FALSE)
    }
    if (is.null(rownames(prob_mat))) rownames(prob_mat) <- cache$valid_cells
    if (is.null(colnames(prob_mat))) colnames(prob_mat) <- cache$valid_cells
  }
  if (!valid_names(rownames(prob_mat), nrow(prob_mat)) ||
      !valid_names(colnames(prob_mat), ncol(prob_mat)) ||
      !setequal(rownames(prob_mat), colnames(prob_mat))) {
    stop("P_norm axes must contain the same unique cell names.", call. = FALSE)
  }
  if (!valid_names(rownames(embedding), nrow(embedding))) {
    stop("The reduction must have unique cell names.", call. = FALSE)
  }

  object_cells <- colnames(seurat_obj)
  meta <- seurat_obj[[]]
  labels <- as.character(meta[object_cells, cluster_col])
  seed_cells <- object_cells[!is.na(labels) & labels == metaclone_id]
  if (length(seed_cells) == 0L) {
    stop(sprintf("Metaclone '%s' was not found in '%s'.", metaclone_id, cluster_col), call. = FALSE)
  }
  plot_cells <- intersect(object_cells, rownames(embedding))
  missing_seeds <- setdiff(seed_cells, intersect(plot_cells, rownames(prob_mat)))
  if (length(missing_seeds) > 0L) {
    stop(sprintf("%d seed cell(s) are missing from P_norm or the requested reduction: %s",
                 length(missing_seeds), paste(utils::head(missing_seeds, 5), collapse = ", ")),
         call. = FALSE)
  }
  coordinates <- embedding[plot_cells, dims, drop = FALSE]
  if (any(!is.finite(coordinates))) {
    stop("The selected embedding coordinates contain non-finite values.", call. = FALSE)
  }

  seed_matrix <- prob_mat[seed_cells, , drop = FALSE]
  entries <- if (inherits(seed_matrix, "sparseMatrix") &&
                 "x" %in% methods::slotNames(seed_matrix)) {
    methods::slot(seed_matrix, "x")
  } else {
    as.numeric(seed_matrix)
  }
  if (any(!is.finite(entries)) || any(entries < 0)) {
    stop("The selected P_norm rows contain non-finite or negative probabilities.", call. = FALSE)
  }

  # Preserve the embedding scripts' mean, seed exclusion, and maximum scaling.
  mean_probability <- as.numeric(Matrix::colMeans(seed_matrix))
  names(mean_probability) <- colnames(prob_mat)
  sister_probability <- mean_probability
  sister_probability[seed_cells] <- 0
  probability_max <- max(sister_probability)
  relative_probability <- sister_probability
  if (probability_max > 0) relative_probability <- sister_probability / probability_max

  unscored_cells <- setdiff(plot_cells, colnames(prob_mat))
  if (length(unscored_cells) > 0L) {
    warning(sprintf("%d plotted cell(s) have no cached probability and are shown as background.",
                    length(unscored_cells)), call. = FALSE)
  }
  plot_df <- data.frame(
    Cell_ID = plot_cells,
    Dim1 = coordinates[, 1], Dim2 = coordinates[, 2],
    MeanProbability = unname(mean_probability[plot_cells]),
    SisterProbability = unname(sister_probability[plot_cells]),
    RelativeProbability = unname(relative_probability[plot_cells]),
    IsSeed = plot_cells %in% seed_cells,
    stringsAsFactors = FALSE,
    row.names = plot_cells
  )
  plot_df$IsSister <- !plot_df$IsSeed & !is.na(plot_df$RelativeProbability) &
    plot_df$RelativeProbability > 0
  list(plot_data = plot_df, seed_cells = seed_cells, probability_max = probability_max)
}
