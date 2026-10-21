#' Plot a pCORAL cell or pseudoclone embedding
#'
#' @description
#' Displays existing embedding coordinates with Arial text and no axes.
#' Cell embeddings use filled points; clone embeddings use open circles.
#'
#' @param seurat_obj A Seurat object containing the requested reduction.
#' @param reduction Reduction to plot. Defaults to \code{"CORAL_cell_embed"};
#'   use \code{"CORAL_clone_embed"} for the pseudoclone embedding.
#' @param group_by Optional metadata column for categorical coloring.
#'   \code{NULL} uses the current Seurat identities.
#' @param dims Two embedding dimensions to display.
#' @param colors Optional color vector, preferably named by group.
#' @param point_size Point size.
#' @param alpha Point opacity, between 0 and 1.
#' @param hollow Use open circles. \code{NULL} selects open circles for
#'   \code{"CORAL_clone_embed"} and filled points for other reductions.
#' @param stroke Outline width for open circles.
#' @param base_size,base_family Base font size and font family.
#' @param title Optional plot title.
#'
#' @details
#' Coordinates and grouping labels are aligned by cell name. All cells in
#' the selected reduction are plotted without subsampling. Numeric metadata
#' are treated as categorical labels, and missing labels remain missing and
#' are shown in gray. Named colors must cover all observed non-missing groups.
#'
#' This function reads only the reduction and cell labels. It does not
#' recompute an embedding, access expression layers, modify the Seurat
#' object, print a figure, or save a file.
#'
#' @return A ggplot object. Its \code{data} contains \code{Cell_ID},
#'   \code{Dim1}, \code{Dim2}, and \code{Group}.
#' @examples
#' \dontrun{
#' p_cell <- PlotCORALEmbedding(seurat_obj)
#' print(p_cell)
#'
#' p_clone <- PlotCORALEmbedding(seurat_obj, reduction = "CORAL_clone_embed")
#' print(p_clone)
#'
#' p_metaclone <- PlotCORALEmbedding(
#'   seurat_obj, reduction = "CORAL_clone_embed", group_by = "metaclone"
#' )
#' print(p_metaclone)
#' }
#' @export
PlotCORALEmbedding <- function(
    seurat_obj,
    reduction = "CORAL_cell_embed",
    group_by = NULL,
    dims = c(1, 2),
    colors = NULL,
    point_size = 1.2,
    alpha = 1,
    hollow = NULL,
    stroke = 0.6,
    base_size = 18,
    base_family = "Arial",
    title = NULL
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("'seurat_obj' must be a Seurat object.", call. = FALSE)
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
  embedding <- embedding[, dims, drop = FALSE]
  cells <- rownames(embedding)
  meta <- seurat_obj[[]]
  if (nrow(embedding) == 0L || is.null(cells) || anyNA(cells) ||
      any(!nzchar(cells)) || anyDuplicated(cells) ||
      !all(cells %in% rownames(meta)) || !all(is.finite(embedding))) {
    stop("The reduction must contain finite coordinates and unique cell names matching the metadata.",
         call. = FALSE)
  }

  if (is.null(group_by)) {
    identities <- Seurat::Idents(seurat_obj)
    if (!all(cells %in% names(identities))) {
      stop("Seurat identities could not be matched to the embedding cells.", call. = FALSE)
    }
    groups <- identities[cells]
    legend_title <- "Identity"
  } else {
    if (!is.character(group_by) || length(group_by) != 1L ||
        is.na(group_by) || !group_by %in% colnames(meta)) {
      stop("'group_by' must name an existing metadata column.", call. = FALSE)
    }
    groups <- meta[cells, group_by]
    legend_title <- group_by
  }
  if (!is.atomic(groups) || !is.null(dim(groups))) {
    stop("Grouping labels must be a one-dimensional atomic vector.", call. = FALSE)
  }
  if (!is.factor(groups)) groups <- factor(groups)
  plot_df <- data.frame(
    Cell_ID = cells, Dim1 = embedding[, 1], Dim2 = embedding[, 2],
    Group = groups, stringsAsFactors = FALSE, row.names = cells
  )

  positive_parameters <- list(point_size = point_size, base_size = base_size)
  for (parameter in names(positive_parameters)) {
    value <- positive_parameters[[parameter]]
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop(sprintf("'%s' must be one positive number.", parameter), call. = FALSE)
    }
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be between 0 and 1.", call. = FALSE)
  }
  if (is.null(hollow)) hollow <- identical(reduction, "CORAL_clone_embed")
  if (!is.logical(hollow) || length(hollow) != 1L || is.na(hollow)) {
    stop("'hollow' must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  if (!is.numeric(stroke) || length(stroke) != 1L ||
      !is.finite(stroke) || stroke < 0) {
    stop("'stroke' must be one non-negative number.", call. = FALSE)
  }
  if (!is.character(base_family) || length(base_family) != 1L || is.na(base_family)) {
    stop("'base_family' must be one font-family name.", call. = FALSE)
  }
  if (!is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("'title' must be NULL or one character string.", call. = FALSE)
  }

  p <- ggplot2::ggplot(
    plot_df, ggplot2::aes(x = .data$Dim1, y = .data$Dim2, color = .data$Group)
  ) +
    ggplot2::geom_point(
      size = point_size, alpha = alpha,
      shape = if (hollow) 1 else 16, stroke = stroke
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_family = base_family, base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family, size = base_size),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "plain"),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::labs(color = legend_title, title = title)

  if (is.null(colors)) {
    p <- p + ggplot2::scale_color_discrete(na.value = "#A9A9A9")
  } else {
    if (!is.character(colors) || length(colors) == 0L || anyNA(colors)) {
      stop("'colors' must be a non-empty character vector of colors.", call. = FALSE)
    }
    observed_groups <- levels(droplevels(plot_df$Group))
    if (is.null(names(colors))) {
      if (length(colors) < length(observed_groups)) {
        stop("Provide at least one color per observed group.", call. = FALSE)
      }
    } else if (anyNA(names(colors)) || any(!nzchar(names(colors))) ||
               anyDuplicated(names(colors)) ||
               !all(observed_groups %in% names(colors))) {
      stop("Named colors must uniquely cover every observed group.", call. = FALSE)
    }
    grDevices::col2rgb(colors)
    p <- p + ggplot2::scale_color_manual(values = colors, na.value = "#A9A9A9")
  }
  p
}
