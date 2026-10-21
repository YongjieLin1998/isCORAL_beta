# ==============================================================================
# Metaclone Analysis & Visualization API Suite (Ultra-Integrated Version)
# ==============================================================================
# Dependency Libraries: ggplot2, dplyr, patchwork, Matrix, ggrepel

# Register data-frame columns used through ggplot2 non-standard evaluation.
utils::globalVariables(c(
  "Dim1", "Dim2", "Cluster", "TargetCluster", "RenderGroup", "Prob",
  "IsSeed", "IsSibling", "U1", "U2", "P1", "P2", "Seed_ID",
  "C1", "C2", "Target_Cluster", "Seed_Lineage", "x", "y", "xend", "yend",
  "CellType", "Metaclone", "Freq", "ExpectedPct", "CapturePct"
))

.resolve_metaclone_algorithm <- function(
    algorithm, verbose = TRUE,
    leidenbase_available = requireNamespace("leidenbase", quietly = TRUE),
    seurat_version = utils::packageVersion("Seurat")
) {
  algorithm_requested <- algorithm
  clustering_engine <- if (algorithm == 1) "Louvain" else sprintf("Seurat algorithm %s", algorithm)

  if (algorithm == 4) {
    seurat_version_text <- as.character(seurat_version)
    seurat_supports_leidenbase <- tryCatch(
      numeric_version(seurat_version_text) >= numeric_version("5.2.0"),
      error = function(e) FALSE
    )

    if (!isTRUE(leidenbase_available)) {
      warning(
        paste0(
          "Package 'leidenbase' is not available; Leiden clustering ",
          "(algorithm = 4) cannot be used. Falling back to Louvain ",
          "(algorithm = 1). isCORAL does not install optional packages ",
          "at runtime; install 'leidenbase' before rerunning to use Leiden."
        ),
        call. = FALSE
      )
      algorithm <- 1
      clustering_engine <- "Louvain fallback (leidenbase unavailable)"
    } else if (!isTRUE(seurat_supports_leidenbase)) {
      warning(
        sprintf(
          paste0(
            "Seurat %s predates the leidenbase-backed Leiden implementation ",
            "available in Seurat >= 5.2.0. Falling back to Louvain ",
            "(algorithm = 1); isCORAL does not invoke the legacy ",
            "Python-backed 'leiden' package."
          ),
          seurat_version_text
        ),
        call. = FALSE
      )
      algorithm <- 1
      clustering_engine <- "Louvain fallback (Seurat < 5.2.0)"
    } else {
      clustering_engine <- "Leiden (leidenbase)"
      if (verbose) message("   Using the C++ Leiden engine from 'leidenbase'.")
    }
  }

  list(
    algorithm_requested = algorithm_requested,
    algorithm_effective = algorithm,
    clustering_engine = clustering_engine
  )
}

#' @title Get a dynamically adaptive topology-aware Metaclone color palette
#' @description Uses high-contrast preset colors for N <= 18 and interpolates the palette for N > 18.
#' @param n Integer number of Metaclone clusters in the current manifold.
#' @param cluster_names Character vector containing the cluster labels used to name the palette.
#' @return A hexadecimal color vector named with the supplied cluster labels.
#' @import ggplot2
#' @import patchwork
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom grid unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats dist hclust
#' @export
GetMetacloneColors <- function(n, cluster_names) {


  base_colors <- c(
    "#FB9A99", "#56B4E9", "#33A02C", "#CAB2D6",
    "#E31A1C", "#B15928", "#FDBF6F", "#E6AB02",
    "#08519C", "#FF7F00", "#6A3D9A", "#E7298A",
    "#B3DE69", "#8DD3C7", "#6B8E23", "#B2DF8A",
    "#A6761D", "#BC80BD"
  )
  if (n <= length(base_colors)) {
    res_colors <- base_colors[1:n]
  } else {
    message(sprintf("   The number of Metaclones (%d) exceeds 18; extending the palette by color interpolation...", n))
    pal_func <- colorRampPalette(base_colors)
    res_colors <- pal_func(n)
  }
  names(res_colors) <- cluster_names
  return(res_colors)
}


#' @title Identify Metaclone communities in the CORAL clone space
#' @description Builds a shared-nearest-neighbor graph in the CORAL clone
#'   embedding and identifies communities with Leiden or Louvain clustering.
#'   Leiden uses the optional \code{leidenbase} package through Seurat >=
#'   5.2.0. If that backend is unavailable, a requested Leiden run warns and
#'   falls back to Louvain.
#' @param seurat_obj A Seurat object containing the requested clone-space reduction.
#' @param resolution Community-detection resolution passed to
#'   \code{Seurat::FindClusters()}.
#' @param reduction Name of the CORAL clone-space reduction.
#' @param dims Dimensions used to construct the nearest-neighbor graph. The
#'   historical default \code{1:30} is limited to the available dimensions for
#'   compatibility with lower-dimensional CORAL reductions; other out-of-range
#'   requests remain errors.
#' @param graph_prefix Prefix for generated nearest-neighbor graph names.
#' @param cluster_col_name Metadata column receiving the Metaclone labels.
#' @param algorithm Clustering algorithm passed to
#'   \code{Seurat::FindClusters()}; 4 requests Leiden and 1 requests Louvain.
#' @param random_seed Random seed for community detection.
#' @param verbose Whether to report graph construction and clustering progress.
#'   Defaults to \code{FALSE}; warnings and errors remain visible.
#' @return The modified Seurat object with Metaclone labels, neighbor graphs,
#'   and a run record in \code{@misc$isCORAL_metaclone}.
#' @details \code{FindMetaclone()} is the canonical beta 1.0.0 entry point.
#'   The deprecated \code{FindMetaclones()} name is retained only as a
#'   compatibility wrapper for scripts written for versions 6.1.0 and 6.1.1.
#'   Optional software is never installed while this function is running.
#' @export
FindMetaclone <- function(
    seurat_obj, resolution = 0.7, reduction = "CORAL_clone_embed",
    dims = 1:30, graph_prefix = "pcoral", cluster_col_name = "metaclone", 
    algorithm = 4, random_seed = 42, verbose = FALSE
) {
  metaclone_call <- match.call()
  
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction '%s' was not found. Run the embedding workflow first.", reduction))
  }
  available_dims <- ncol(Seurat::Embeddings(seurat_obj, reduction))
  if (!is.numeric(dims) || length(dims) == 0 || anyNA(dims) ||
      any(!is.finite(dims)) || any(dims < 1) || any(dims != floor(dims))) {
    stop("'dims' must be a non-empty vector of positive integer indices.")
  }
  dims <- as.integer(dims)
  requested_dims <- dims
  if (any(dims > available_dims)) {
    is_historical_default <- length(dims) == 30L &&
      identical(dims, seq_len(30L))
    if (!is_historical_default) {
      stop(sprintf(
        "The requested dimensions exceed reduction '%s', which contains %d dimensions.",
        reduction, available_dims
      ))
    }
    dims <- dims[dims <= available_dims]
    if (length(dims) == 0L) {
      stop(sprintf("Reduction '%s' does not contain a usable dimension.", reduction))
    }
    warning(sprintf(
      "The historical default dims = 1:30 was limited to 1:%d for reduction '%s'.",
      available_dims, reduction
    ), call. = FALSE)
  }
  
  algorithm_resolution <- .resolve_metaclone_algorithm(
    algorithm = algorithm,
    verbose = verbose
  )
  algorithm_requested <- algorithm_resolution$algorithm_requested
  algorithm <- algorithm_resolution$algorithm_effective
  clustering_engine <- algorithm_resolution$clustering_engine
  
  nn_name  <- paste0(graph_prefix, "_nn")
  snn_name <- paste0(graph_prefix, "_snn")
  
  if (verbose) message(sprintf(
    ">>> [1/2] Building the topology graph (Reduction: %s, Dims: %s)...",
    reduction, paste(dims, collapse = ",")
  ))
  seurat_obj <- .isCORAL_with_seurat_output(
    Seurat::FindNeighbors(
      seurat_obj, reduction = reduction, dims = dims,
      graph.name = c(nn_name, snn_name), verbose = FALSE
    ), verbose = verbose
  )
  
  alg_name <- ifelse(algorithm == 4, "Leiden", ifelse(algorithm == 1, "Louvain", "Other"))
  if (verbose) message(sprintf(">>> [2/2] Identifying Metaclone communities (Algorithm: %s, Resolution: %s)...", alg_name, resolution))
  
  find_clusters_args <- list(
    object = seurat_obj,
    graph.name = snn_name,
    resolution = resolution,
    algorithm = algorithm,
    random.seed = random_seed,
    verbose = verbose
  )
  if (algorithm == 4) find_clusters_args$leiden_method <- "leidenbase"
  seurat_obj <- .isCORAL_with_seurat_output(
    do.call(Seurat::FindClusters, find_clusters_args), verbose = verbose
  )
  
  active_idents <- Seurat::Idents(seurat_obj)
  cluster_labels <- as.character(active_idents)
  if (length(cluster_labels) != ncol(seurat_obj)) {
    stop("The clustering result does not contain exactly one label for every cell.")
  }
  if (is.null(names(cluster_labels))) names(cluster_labels) <- names(active_idents)
  if (is.null(names(cluster_labels))) names(cluster_labels) <- colnames(seurat_obj)
  if (anyNA(names(cluster_labels)) || any(!nzchar(names(cluster_labels))) ||
      anyDuplicated(names(cluster_labels)) ||
      !setequal(colnames(seurat_obj), names(cluster_labels))) {
    stop("The clustering result could not be aligned to all cells by cell barcode.")
  }
  seurat_obj[[cluster_col_name]] <- cluster_labels[colnames(seurat_obj)]
  Seurat::Idents(seurat_obj) <- cluster_col_name
  
  package_version_or_na <- function(package_name) {
    tryCatch(
      as.character(utils::packageVersion(package_name)),
      error = function(e) NA_character_
    )
  }
  run_timestamp_utc <- format(
    as.POSIXct(Sys.time(), tz = "UTC"),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )
  cluster_values <- as.character(seurat_obj[[cluster_col_name]][[1]])
  cluster_counts <- table(cluster_values, useNA = "ifany")
  run_record <- list(
    entry_point = "FindMetaclone",
    call = metaclone_call,
    call_text = paste(deparse(metaclone_call), collapse = " "),
    timestamp_utc = run_timestamp_utc,
    requested_parameters = list(
      resolution = resolution,
      reduction = reduction,
      dims = requested_dims,
      graph_prefix = graph_prefix,
      cluster_col_name = cluster_col_name,
      algorithm = algorithm_requested,
      random_seed = random_seed,
      verbose = verbose
    ),
    effective_parameters = list(
      resolution = resolution,
      reduction = reduction,
      dims = dims,
      available_dimensions = available_dims,
      graph_prefix = graph_prefix,
      cluster_col_name = cluster_col_name,
      algorithm = algorithm,
      clustering_engine = clustering_engine,
      leiden_method = if (algorithm == 4) "leidenbase" else NULL,
      random_seed = random_seed,
      verbose = verbose
    ),
    outputs = list(
      nearest_neighbor_graph = nn_name,
      shared_nearest_neighbor_graph = snn_name,
      cluster_metadata_column = cluster_col_name,
      active_identity = cluster_col_name,
      cluster_counts = cluster_counts,
      number_of_clusters = length(cluster_counts),
      number_of_cells = ncol(seurat_obj)
    ),
    provenance = list(
      iscoral_version = "1.0.0",
      release_status = "beta",
      R = R.version.string,
      Seurat = package_version_or_na("Seurat"),
      SeuratObject = package_version_or_na("SeuratObject"),
      leidenbase = package_version_or_na("leidenbase")
    )
  )

  previous_cache <- seurat_obj@misc$isCORAL_metaclone
  if (is.null(previous_cache)) {
    run_history <- list(run_record)
  } else if (is.list(previous_cache) && !is.null(previous_cache$history)) {
    run_history <- c(previous_cache$history, list(run_record))
  } else {
    run_history <- list(previous_cache, run_record)
  }
  seurat_obj@misc$isCORAL_metaclone <- list(
    schema_version = "1.0",
    latest = run_record,
    history = run_history
  )

  return(seurat_obj)
}

#' @rdname FindMetaclone
#' @export
FindMetaclones <- function(
    seurat_obj, resolution = 0.7, reduction = "CORAL_clone_embed",
    dims = 1:30, graph_prefix = "pcoral", cluster_col_name = "metaclone",
    algorithm = 4, random_seed = 42, verbose = FALSE
) {
  legacy_call <- match.call()
  .Deprecated("FindMetaclone", package = "isCORAL")
  result <- FindMetaclone(
    seurat_obj = seurat_obj,
    resolution = resolution,
    reduction = reduction,
    dims = dims,
    graph_prefix = graph_prefix,
    cluster_col_name = cluster_col_name,
    algorithm = algorithm,
    random_seed = random_seed,
    verbose = verbose
  )
  result@misc$isCORAL_metaclone$latest$entry_point <- "FindMetaclones"
  result@misc$isCORAL_metaclone$latest$compatibility_target <- "FindMetaclone"
  result@misc$isCORAL_metaclone$latest$call <- legacy_call
  result@misc$isCORAL_metaclone$latest$call_text <-
    paste(deparse(legacy_call), collapse = " ")
  history_index <- length(result@misc$isCORAL_metaclone$history)
  result@misc$isCORAL_metaclone$history[[history_index]] <-
    result@misc$isCORAL_metaclone$latest
  result
}

#' @title Plot the Metaclone landscape and faceted density contours
#' @description Draws a landscape panel and horizontal cluster facets using a consistent open-circle style.
#' @param seurat_obj A Seurat object containing the requested reduction and Metaclone labels.
#' @param reduction Reduction used for the two-dimensional landscape.
#' @param cluster_col Metadata column containing Metaclone labels.
#' @param custom_colors Optional named vector of cluster colors.
#' @return A patchwork object containing the landscape and cluster facets.
#' @export
PlotMetacloneLandscape <- function(seurat_obj, reduction = "CORAL_clone_embed", cluster_col = "metaclone", custom_colors = NULL) {
  mat <- Seurat::Embeddings(seurat_obj, reduction)
  df <- data.frame(
    Dim1 = mat[, 1], Dim2 = mat[, 2],
    Cluster = as.character(seurat_obj[[cluster_col]][[1]]),
    stringsAsFactors = FALSE
  )
  
  cluster_levels <- as.character(sort(as.numeric(unique(df$Cluster))))
  df$Cluster <- factor(df$Cluster, levels = cluster_levels)
  
  if (is.null(custom_colors)) {
    custom_colors <- GetMetacloneColors(n = length(cluster_levels), cluster_names = cluster_levels)
  }
  
  # --- Upper panel: landscape ---
  p_top <- ggplot(df, aes(x = Dim1, y = Dim2, color = Cluster)) +
    # Use open circles and a controlled stroke width for a consistent appearance.
    geom_point(alpha = 0.7, size = 1.2, shape = 1, stroke = 0.6) + 
    geom_density_2d(linewidth = 1.0, alpha = 1.0) +
    scale_color_manual(values = custom_colors) + 
    theme_classic(base_family = "Arial") +
    theme(legend.position = "right", legend.title = element_text(face = "plain"),
          axis.title = element_text(size = 14), axis.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "plain", margin = margin(b=10))) +
    labs(x = paste0(reduction, "_1"), y = paste0(reduction, "_2"), color = "clone community", title = "metaclone embedding")
  
  # --- Lower panel: cluster facets ---
  df_bg <- df[, c("Dim1", "Dim2")]
  p_bottom <- ggplot() +
    # Draw background cells as fine open circles.
    geom_point(data = df_bg, aes(x = Dim1, y = Dim2), color = "#E0E4E7", size = 0.3, shape = 1, stroke = 0.3) +
    # Use the same open-circle style for highlighted cells.
    geom_point(data = df, aes(x = Dim1, y = Dim2, color = Cluster), alpha = 0.6, size = 0.6, shape = 1, stroke = 0.5) +
    geom_density_2d(data = df, aes(x = Dim1, y = Dim2, color = Cluster), linewidth = 0.8, alpha = 1.0) +
    scale_color_manual(values = custom_colors) + 
    facet_wrap(~ Cluster, nrow = 2) + 
    theme_classic(base_family = "Arial") +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey95", color = "grey80"),
          strip.text = element_text(size = 11, face = "plain"), axis.text = element_blank(), 
          axis.ticks = element_blank(), axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "plain", margin = margin(b=10))) +
    labs(title = "metaclone separated")
  
  p_combined <- p_top / p_bottom + 
    plot_layout(heights = c(1.5, 1)) +
    plot_annotation(
      title = sprintf("metaclone clustering"),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "plain", margin = margin(b=15)))
    )
  return(p_combined)
}

#' @title Plot a validation matrix of inferred and experimental clone structure
#' @description Draws a complete faceted validation matrix with seed cells shown in darker cluster colors.
#' @param seurat_obj A Seurat object containing clone-space coordinates,
#'   Metaclone labels, and cached pseudoclone probabilities.
#' @param reduction Reduction used for the validation panels.
#' @param space_title Text identifying the plotted space.
#' @param barcode_col Required metadata column containing experimental ground-truth clone labels.
#' @param cluster_col Metadata column containing Metaclone labels.
#' @param na_barcodes Optional character vector of additional barcode labels treated as invalid or background.
#' @param custom_colors Optional named vector of cluster colors.
#' @return A patchwork object comparing inferred pseudoclone structure with
#'   experimental clone structure.
#' @export
PlotValidationSuite <- function(
    seurat_obj, reduction, space_title, 
    barcode_col, cluster_col = "metaclone", 
    na_barcodes = NULL,   # Optional additions to the default invalid-barcode set.
    custom_colors = NULL
) {
  # Require rasterized point layers without modifying the user's library.
  if (!requireNamespace("ggrastr", quietly = TRUE)) stop("Package 'ggrastr' is required for PlotValidationSuite().")
  
  if (missing(barcode_col)) stop("barcode_col is required and must identify the experimental clone-label column.")
  if (is.null(seurat_obj@misc$isCORAL_embedding$P_norm)) stop("P_norm was not found in seurat_obj@misc$isCORAL_embedding.")
  
  # ==============================================================================
  # Combine user-supplied invalid labels with the default set.
  # ==============================================================================
  # "0" is included as a default invalid barcode.
  base_na_barcodes <- c("NA", "Unknown_NA", NA, "", "0")
  final_na_barcodes <- unique(c(base_na_barcodes, na_barcodes))
  
  mat <- Seurat::Embeddings(seurat_obj, reduction)
  df_base <- data.frame(
    cell = colnames(seurat_obj), Dim1 = mat[, 1], Dim2 = mat[, 2],
    Cluster = as.character(seurat_obj[[cluster_col]][[1]]),
    stringsAsFactors = FALSE
  )
  
  all_energy_clusters <- seurat_obj[[cluster_col]][[1]]
  all_barcodes <- as.character(seurat_obj[[barcode_col]][[1]])
  P_norm <- seurat_obj@misc$isCORAL_embedding$P_norm
  
  cluster_levels <- as.character(sort(as.numeric(unique(all_energy_clusters))))
  if (is.null(custom_colors)) {
    custom_colors <- GetMetacloneColors(n = length(cluster_levels), cluster_names = cluster_levels)
  }
  
  # ==============================================================================
  # Build and combine the base and darkened seed palettes.
  # ==============================================================================
  darken_color <- function(hex, factor = 0.45) {
    rgb_val <- as.numeric(col2rgb(hex)) * factor
    rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255)
  }
  seed_colors <- sapply(custom_colors, darken_color)
  names(seed_colors) <- paste0(names(custom_colors), "_seed") # For example, "1_seed".
  
  # Combine bright cluster colors with darker seed colors.
  full_palette <- c(custom_colors, seed_colors)
  
  get_clean_clones <- function(cluster_id) {
    bcs <- all_barcodes[colnames(seurat_obj) %in% colnames(seurat_obj)[all_energy_clusters == cluster_id]]
    # Filter using the combined invalid-barcode set.
    return(unique(bcs[!(bcs %in% final_na_barcodes) & !is.na(bcs)]))
  }
  
  long_list <- list()
  for (k in sort(unique(all_energy_clusters))) {
    idx <- which(all_energy_clusters == k)
    if (length(idx) == 0) next
    
    prob <- as.numeric(Matrix::colMeans(P_norm[idx, , drop = FALSE]))
    prob[is.na(prob)] <- 0
    if (max(prob) > 0) prob = prob / max(prob)
    
    clones_pool <- get_clean_clones(k)
    siblings <- colnames(seurat_obj)[all_barcodes %in% clones_pool]
    
    tmp <- df_base[, c("cell", "Dim1", "Dim2", "Cluster")]
    tmp$TargetCluster <- as.character(k)
    tmp$Prob <- prob
    tmp$IsSibling <- (tmp$cell %in% siblings)
    tmp$IsSeed <- (tmp$Cluster == as.character(k))
    long_list[[length(long_list) + 1]] <- tmp
  }
  
  df_all <- dplyr::bind_rows(long_list)
  df_all$TargetCluster <- factor(df_all$TargetCluster, levels = cluster_levels)
  
  # ==============================================================================
  # Use a dedicated render group to keep facet colors stable.
  # ==============================================================================
  df_all$RenderGroup <- as.character(df_all$TargetCluster)
  df_all$RenderGroup[df_all$IsSeed] <- paste0(df_all$TargetCluster[df_all$IsSeed], "_seed")
  
  df_inferred <- df_all[order(df_all$Prob, df_all$IsSeed), ]
  df_physical <- df_all[order(df_all$IsSibling, df_all$IsSeed), ]
  
  base_theme <- theme_classic(base_family = "Arial") +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey95", color = "grey80"),
          strip.text = element_text(size = 11, face = "plain"), axis.text = element_blank(), 
          axis.ticks = element_blank(), axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "plain", margin = margin(b=10)))
  
  # ==============================================================================
  # Render the inferred panel with rasterized high-density layers.
  # ==============================================================================
  p_inferred <- ggplot() +
    # Rasterize the high-density low-probability background.
    ggrastr::rasterise(geom_point(data = subset(df_inferred, Prob <= 0.02 & !IsSeed), aes(x = Dim1, y = Dim2), color = "#E0E4E7", size = 0.3), dpi = 300) +
    
    # Rasterize the alpha-mapped probability cloud for efficient PDF rendering.
    ggrastr::rasterise(geom_point(data = subset(df_inferred, Prob > 0.02 & !IsSeed), aes(x = Dim1, y = Dim2, color = RenderGroup, alpha = Prob), size = 0.6), dpi = 300) +
    
    # Keep the small seed-cell layer as vectors.
    geom_point(data = subset(df_inferred, IsSeed), aes(x = Dim1, y = Dim2, color = RenderGroup), shape = 16, size = 1.2, alpha = 1.0) +
    
    scale_color_manual(values = full_palette) + 
    scale_alpha_continuous(range = c(0, 1.0), guide = "none") + 
    facet_wrap(~ TargetCluster, nrow = 2) + base_theme + labs(title = "inferred metaclones")
  
  # ==============================================================================
  # Render the experimental ground-truth panel with rasterized dense layers.
  # ==============================================================================
  p_physical <- ggplot() +
    # Rasterize gray cells outside the target clone.
    ggrastr::rasterise(geom_point(data = subset(df_physical, !IsSibling & !IsSeed), aes(x = Dim1, y = Dim2), color = "#E0E4E7", size = 0.3), dpi = 300) +
    
    # Rasterize cells in the experimental clone distribution.
    ggrastr::rasterise(geom_point(data = subset(df_physical, IsSibling & !IsSeed), aes(x = Dim1, y = Dim2, color = RenderGroup), size = 0.7, alpha = 0.8), dpi = 300) +
    
    # Keep seed cells as vectors.
    geom_point(data = subset(df_physical, IsSeed), aes(x = Dim1, y = Dim2, color = RenderGroup), shape = 16, size = 1.2, alpha = 1.0) +
    
    scale_color_manual(values = full_palette) + 
    facet_wrap(~ TargetCluster, nrow = 2) + base_theme + labs(title = "experimental lineage groundtruth")
  
  p_combined_facets <- p_inferred / p_physical + 
    patchwork::plot_annotation(
      title = sprintf("validation (%s)", space_title),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "plain", margin = margin(b=15)))
    )
  
  return(p_combined_facets)
}

#' @title Plot experimental siblings and inferred probabilities for selected cells
#' @description Contrasts isCORAL pseudoclone probabilities with experimentally
#'   labeled sibling cells for selected seed cells.
#' @param seurat_obj A Seurat object containing the pseudoclone cache and a
#'   \code{umap.raw} reduction.
#' @param demo_cells Character vector of seed-cell names.
#' @param highlight_colors Named colors for the displayed cell lineages.
#' @param clone_col Metadata column containing experimental clone labels.
#' @param bg_gray Color used for background cells.
#' @param point_size Base point size.
#' @return A vertically arranged patchwork object with one
#'   inferred/experimental diptych per valid seed cell.
#' @export
PlotProbabilityCollapseDiptych <- function(
    seurat_obj, demo_cells, highlight_colors, 
    clone_col, bg_gray = "#E0E4E7", point_size = 0.5
) {
  if (missing(clone_col)) stop("clone_col is required and must identify the experimental clone-label column.")
  
  cache <- seurat_obj@misc$isCORAL_embedding
  P_norm <- cache$P_norm
  valid_cells <- cache$valid_cells
  rownames(P_norm) <- valid_cells
  
  mat_umap <- Seurat::Embeddings(seurat_obj, "umap.raw")[valid_cells, 1:2]
  all_clones <- as.character(seurat_obj@meta.data[valid_cells, clone_col])
  names(all_clones) <- valid_cells
  
  theme_collapse <- theme_void(base_family = "Arial") + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "plain"), 
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40", face = "plain"),
      legend.position = "none", 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5) 
    )
  
  plot_list <- list()
  for (i in seq_along(demo_cells)) {
    seed_cell <- demo_cells[i]
    if (!seed_cell %in% valid_cells) next
    
    seed_clone <- all_clones[seed_cell]
    if (is.na(seed_clone)) seed_clone <- "Unknown/NA"
    
    cell_probs <- as.numeric(P_norm[seed_cell, ])
    cell_probs[is.na(cell_probs)] <- 0
    if (max(cell_probs) > 0) cell_probs <- cell_probs / max(cell_probs) 
    
    raw_clusters <- as.character(seurat_obj@meta.data[valid_cells, "cluster_integrated"])
    mapped_clusters <- ifelse(raw_clusters %in% c("HSC", "GMP", "Mk"), raw_clusters, "Other")
    
    df <- data.frame(
      U1 = mat_umap[, 1], U2 = mat_umap[, 2], Prob = cell_probs, Cluster = mapped_clusters,
      IsSeed = (valid_cells == seed_cell),
      IsSibling = (all_clones == seed_clone) & (valid_cells != seed_cell) & (seed_clone != "Unknown/NA"),
      stringsAsFactors = FALSE
    )
    
    df_seed <- df[df$IsSeed, ]
    df_sibling <- df[df$IsSibling, ]
    df_cloud <- df[order(df$Prob), ] 
    
    p_cloud <- ggplot() +
      geom_point(data = subset(df_cloud, Prob == 0), aes(x = U1, y = U2), color = bg_gray, size = point_size) +
      geom_point(data = subset(df_cloud, Prob > 0), aes(x = U1, y = U2, color = Cluster, alpha = Prob), size = point_size * 1.6) +
      geom_point(data = df_seed, aes(x = U1, y = U2, fill = Cluster), shape = 21, size = 5.0, color = "black", stroke = 1.2) +
      scale_color_manual(values = highlight_colors) + scale_fill_manual(values = highlight_colors) +
      scale_alpha_continuous(range = c(0, 1), guide = "none") + theme_collapse +
      labs(title = "isCORAL inference", subtitle = sprintf("seed cell type: %s", df_seed$Cluster))
    
    clean_clone_id <- sub("^[Ss]apphire_?", "", seed_clone)
    p_gt <- ggplot() +
      geom_point(data = subset(df, !IsSibling & !IsSeed), aes(x = U1, y = U2), color = bg_gray, size = point_size, alpha = 0.8) +
      geom_point(data = df_sibling, aes(x = U1, y = U2, fill = Cluster), shape = 21, size = 2.8, color = "black", stroke = 0.7) +
      geom_point(data = df_seed, aes(x = U1, y = U2, fill = Cluster), shape = 21, size = 5.0, color = "black", stroke = 1.2) +
      scale_fill_manual(values = highlight_colors) + theme_collapse +
      labs(title = "LARRY groundtruth", subtitle = sprintf("Clone ID: %s", clean_clone_id))
    
    plot_list[[i]] <- (p_cloud | p_gt)
  }
  return(wrap_plots(plot_list, ncol = 1))
}


#' @title Plot global fate-blended anchors in the clone space
#' @description Blends lineage colors according to inferred pseudoclone
#'   probabilities for selected seed cells.
#' @param seurat_obj A Seurat object containing the requested reduction and
#'   pseudoclone cache.
#' @param demo_cells Character vector of seed-cell names.
#' @param highlight_colors Named colors for HSC, GMP, Mk, and Other lineages.
#' @param reduction Reduction used to display clone-space coordinates.
#' @return A ggplot object.
#' @export
PlotFateBlendedAnchorMap <- function(seurat_obj, demo_cells, highlight_colors, reduction = "CORAL_clone_embed") {
  cache_global <- seurat_obj@misc$isCORAL_embedding
  valid_cells_global <- cache_global$valid_cells
  
  hex_to_rgb <- function(hex) { as.numeric(col2rgb(hex)) / 255 }
  rgb_HSC   <- hex_to_rgb(highlight_colors["HSC"])
  rgb_GMP   <- hex_to_rgb(highlight_colors["GMP"])
  rgb_Mk    <- hex_to_rgb(highlight_colors["Mk"])
  rgb_Other <- hex_to_rgb(highlight_colors["Other"]) 
  
  mat_clone_global <- Seurat::Embeddings(seurat_obj, reduction)[valid_cells_global, 1:2]
  raw_clusters_global <- as.character(seurat_obj@meta.data[valid_cells_global, "cluster_integrated"])
  mapped_clusters_global <- ifelse(raw_clusters_global %in% c("HSC", "GMP", "Mk"), raw_clusters_global, "Other")
  
  P_norm_global <- cache_global$P_norm
  rownames(P_norm_global) <- valid_cells_global
  
  is_HSC   <- mapped_clusters_global == "HSC"
  is_GMP   <- mapped_clusters_global == "GMP"
  is_Mk    <- mapped_clusters_global == "Mk"
  is_Other <- mapped_clusters_global == "Other"
  
  prob_HSC_all   <- as.numeric(Matrix::rowSums(P_norm_global[, is_HSC, drop = FALSE]))
  prob_GMP_all   <- as.numeric(Matrix::rowSums(P_norm_global[, is_GMP, drop = FALSE]))
  prob_Mk_all    <- as.numeric(Matrix::rowSums(P_norm_global[, is_Mk, drop = FALSE]))
  prob_Other_all <- as.numeric(Matrix::rowSums(P_norm_global[, is_Other, drop = FALSE]))
  
  total_core_prob <- prob_HSC_all + prob_GMP_all + prob_Mk_all + prob_Other_all
  total_core_prob[total_core_prob == 0] <- 1e-10 
  
  w_HSC   <- prob_HSC_all   / total_core_prob
  w_GMP   <- prob_GMP_all   / total_core_prob
  w_Mk    <- prob_Mk_all    / total_core_prob
  w_Other <- prob_Other_all / total_core_prob
  
  df_global <- data.frame(
    P1 = mat_clone_global[, 1], P2 = mat_clone_global[, 2],
    Cluster = mapped_clusters_global, stringsAsFactors = FALSE
  )
  rownames(df_global) <- valid_cells_global
  
  df_global$Blended_Color <- rgb(
    red   = pmin(pmax(w_HSC * rgb_HSC[1] + w_GMP * rgb_GMP[1] + w_Mk * rgb_Mk[1] + w_Other * rgb_Other[1], 0), 1),
    green = pmin(pmax(w_HSC * rgb_HSC[2] + w_GMP * rgb_GMP[2] + w_Mk * rgb_Mk[2] + w_Other * rgb_Other[2], 0), 1),
    blue  = pmin(pmax(w_HSC * rgb_HSC[3] + w_GMP * rgb_GMP[3] + w_Mk * rgb_Mk[3] + w_Other * rgb_Other[3], 0), 1)
  )
  
  df_seeds_combined <- df_global[demo_cells, ]
  df_seeds_combined$Seed_ID <- paste0("Seed ", 1:nrow(df_seeds_combined))
  
  p <- ggplot() +
    geom_point(data = df_global, aes(x = P1, y = P2), 
               shape = 1, color = df_global$Blended_Color, size = 1.2, stroke = 0.5, alpha = 0.85) +
    geom_point(data = df_seeds_combined, aes(x = P1, y = P2, color = Cluster), 
               shape = 21, fill = "white", size = 6.5, stroke = 2.5) +
    geom_text_repel(data = df_seeds_combined, aes(x = P1, y = P2, label = Seed_ID),
                    size = 5, fontface = "plain", box.padding = 0.8, 
                    point.padding = 0.5, segment.color = "grey50") +
    scale_color_manual(values = highlight_colors) + theme_void(base_family = "Arial") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "plain"), 
          legend.position = "right", legend.title = element_text(size = 12, face = "plain"), 
          legend.text = element_text(size = 11, face = "plain")) +
    labs(title = "CORAL clone embedding", subtitle = "Fate blended pseudoclone space", color = "seed cell lineage")
  return(p)
}


#' @title Plot faceted continuous probability gradients with blended colors
#' @description Displays normalized pseudoclone probabilities for selected seed
#'   cells across the cached cell embedding.
#' @param seurat_obj A Seurat object containing the pseudoclone cache.
#' @param demo_cells Character vector of seed-cell names.
#' @param highlight_colors Named colors for the displayed target lineages.
#' @param bg_gray Color used for zero-probability cells.
#' @return A faceted ggplot object.
#' @export
PlotProbabilityGradientFacets <- function(seurat_obj, demo_cells, highlight_colors, bg_gray = "#E0E4E7") {
  cache_global <- seurat_obj@misc$isCORAL_embedding
  valid_cells_global <- cache_global$valid_cells
  mat_cell <- Seurat::Embeddings(seurat_obj, "CORAL_cell_embed")[valid_cells_global, 1:2]
  raw_clusters_global <- as.character(seurat_obj@meta.data[valid_cells_global, "cluster_integrated"])
  mapped_clusters_global <- ifelse(raw_clusters_global %in% c("HSC", "GMP", "Mk"), raw_clusters_global, "Other")
  
  seed_lineages <- c("HSC", "GMP", "Mk")
  df_facet_combined <- data.frame()
  
  rownames(cache_global$P_norm) <- valid_cells_global
  
  for (i in seq_along(demo_cells)) {
    seed_cell <- demo_cells[i]
    lineage <- seed_lineages[i]
    
    probs <- as.numeric(cache_global$P_norm[seed_cell, ])
    probs[is.na(probs)] <- 0
    if (max(probs) > 0) probs <- probs / max(probs) 
    
    temp_df <- data.frame(
      C1 = mat_cell[, 1], C2 = mat_cell[, 2], Prob = probs,
      Seed_Lineage = factor(lineage, levels = c("HSC", "GMP", "Mk")), 
      Target_Cluster = mapped_clusters_global,                        
      IsSeed = (valid_cells_global == seed_cell)
    )
    temp_df <- temp_df[order(temp_df$Prob), ]
    df_facet_combined <- rbind(df_facet_combined, temp_df)
  }
  
  p <- ggplot(df_facet_combined, aes(x = C1, y = C2)) +
    geom_point(data = subset(df_facet_combined, Prob == 0), color = bg_gray, size = 0.5) +
    geom_point(data = subset(df_facet_combined, Prob > 0), aes(color = Target_Cluster, alpha = Prob), size = 0.8) +
    geom_point(data = subset(df_facet_combined, IsSeed == TRUE), aes(fill = Target_Cluster), shape = 21, color = "black", size = 5.0, stroke = 1.2) +
    scale_color_manual(values = highlight_colors, name = "target cell lineage") +
    scale_fill_manual(values = highlight_colors, guide = "none") +
    scale_alpha_continuous(range = c(0, 1), name = "pseudoclone probability") +
    facet_wrap(~ Seed_Lineage, ncol = 3) + theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "plain", margin = margin(b = 15)),
      strip.text = element_text(size = 14, face = "plain", margin = margin(b = 10)),
      legend.position = "bottom",
      legend.title = element_text(size = 15, face = "plain", vjust = 1),
      legend.text = element_text(size = 14, face = "plain"),
      legend.key.width = unit(3, "line"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.spacing = unit(1.5, "lines"),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    ) +
    labs(title = "Pseudoclone probability distributions on cell embedding")
  return(p)
}

# ==============================================================================
# Dendrogram-enhanced symmetric matrix heatmap API suite
# ==============================================================================

#' @title Plot the experimental seed heatmap with composition clustering
#' @description Summarizes the observed seed-cell composition of each Metaclone.
#' @param seurat_obj A Seurat object containing Metaclone and cell-type metadata.
#' @param cluster_col Metadata column containing Metaclone labels.
#' @param cell_type_col Metadata column containing cell-type labels.
#' @param target_lineages Optional ordered character vector of lineages to display.
#' @return A list with the clustered count table and combined plot.
#' @export
PlotMetacloneSeedHeatmap <- function(
    seurat_obj, cluster_col = "metaclone", cell_type_col = "cluster_integrated", target_lineages = NULL
) {
  if (!requireNamespace("ggdendro", quietly = TRUE)) stop("Package 'ggdendro' is required for PlotMetacloneSeedHeatmap().")
  
  all_clusters <- as.character(seurat_obj[[cluster_col]][[1]])
  all_cell_types <- as.character(seurat_obj[[cell_type_col]][[1]])
  unique_clusters <- as.character(sort(as.numeric(unique(all_clusters))))
  if (is.null(target_lineages)) target_lineages <- sort(unique(all_cell_types))
  
  seed_table <- table(Metaclone = all_clusters, CellType = all_cell_types)
  seed_mat <- as.matrix(seed_table[unique_clusters, target_lineages, drop = FALSE])
  seed_mat[is.na(seed_mat)] <- 0
  
  # Hierarchical clustering
  seed_prop <- prop.table(seed_mat, margin = 1)
  seed_prop[is.na(seed_prop)] <- 0
  hc <- hclust(dist(seed_prop))
  clustered_order <- rownames(seed_mat)[hc$order]
  seed_mat <- seed_mat[clustered_order, , drop = FALSE]
  
  # Console output
  cat("\n========================================================================\n")
  cat("📊 [TRINITY 1/3]: seed cell type\n")
  cat("========================================================================\n")
  print(seed_mat)
  cat("========================================================================\n\n")
  
  # Color scale
  metaclone_colors <- GetMetacloneColors(n = length(unique_clusters), cluster_names = unique_clusters)
  axis_colors <- metaclone_colors[rev(clustered_order)]
  
  # Extract dendrogram data
  dhc <- ggdendro::dendro_data(hc, type = "rectangle")
  N <- length(unique_clusters)
  
  p_dendro <- ggplot(ggdendro::segment(dhc)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5, color = "grey40") +
    coord_flip() + scale_x_continuous(limits = c(N + 0.5, 0.5), expand = c(0, 0)) +
    scale_y_reverse(expand = expansion(mult = c(0.1, 0.02))) + theme_void() + theme(plot.margin = margin(0, -5, 0, 0))
  
  # Build the heatmap
  df_seed <- reshape2::melt(seed_mat)
  colnames(df_seed) <- c("Metaclone", "CellType", "Freq")
  df_seed$Metaclone <- factor(df_seed$Metaclone, levels = rev(clustered_order)) 
  df_seed$CellType <- factor(df_seed$CellType, levels = target_lineages)
  
  p_main <- ggplot(df_seed, aes(x = CellType, y = Metaclone, fill = Freq)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(Freq > 0, Freq, "")), size = 3.2, color = "black", family = "Arial") +
    scale_fill_gradient(low = "#F4F7F9", high = "#1E88E5", name = "Seed Cells") +
    scale_y_discrete(expand = c(0, 0.5)) +
    theme_minimal(base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
      axis.text.y = element_text(size = 12, color = axis_colors, face = "bold"),
      axis.title = element_text(size = 13), panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 15, face = "plain", margin = margin(b=12)), plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(x = "Origin Cell Identity", y = "Metaclone ID", title = "seed cell type")
  
  p_combined <- p_dendro + p_main + plot_layout(widths = c(1, 8))
  return(list(table = seed_mat, plot = p_combined))
}


#' @title Plot the inferred Markov-potential heatmap with a clustering tree
#' @description Aggregates cached pseudoclone probabilities by Metaclone and
#'   target lineage.
#' @param seurat_obj A Seurat object containing Metaclone labels, cell types,
#'   and cached pseudoclone probabilities.
#' @param cluster_col Metadata column containing Metaclone labels.
#' @param cell_type_col Metadata column containing cell-type labels.
#' @param target_lineages Optional ordered character vector of lineages to display.
#' @return A list with the inferred percentage table and combined plot.
#' @export
PlotMetacloneInferredHeatmap <- function(
    seurat_obj, cluster_col = "metaclone", cell_type_col = "cluster_integrated", target_lineages = NULL
) {
  if (is.null(seurat_obj@misc$isCORAL_embedding$P_norm)) stop("P_norm was not found in seurat_obj@misc$isCORAL_embedding.")
  
  all_clusters <- as.character(seurat_obj[[cluster_col]][[1]])
  all_cell_types <- as.character(seurat_obj[[cell_type_col]][[1]])
  P_norm <- seurat_obj@misc$isCORAL_embedding$P_norm
  unique_clusters <- as.character(sort(as.numeric(unique(all_clusters))))
  if (is.null(target_lineages)) target_lineages <- sort(unique(all_cell_types))
  
  inferred_matrix <- matrix(0, nrow = length(unique_clusters), ncol = length(target_lineages))
  rownames(inferred_matrix) <- unique_clusters
  colnames(inferred_matrix) <- target_lineages
  
  for (k in unique_clusters) {
    idx <- which(all_clusters == k)
    if (length(idx) == 0) next
    
    prob_vec <- as.numeric(Matrix::colMeans(P_norm[idx, , drop = FALSE]))
    prob_vec[is.na(prob_vec)] <- 0
    prob_vec[idx] <- 0  # Exclude the seed cluster itself.
    
    raw_sums <- tapply(prob_vec, all_cell_types, sum)
    raw_sums[is.na(raw_sums)] <- 0
    
    total_prob <- sum(raw_sums)
    rel_sums <- if (total_prob > 0) (raw_sums / total_prob) * 100 else raw_sums
    for (lin in target_lineages) if (lin %in% names(rel_sums)) inferred_matrix[k, lin] <- round(as.numeric(rel_sums[lin]), 1)
  }
  
  # Hierarchical clustering
  hc <- hclust(dist(inferred_matrix))
  clustered_order <- rownames(inferred_matrix)[hc$order]
  inferred_matrix <- inferred_matrix[clustered_order, , drop = FALSE]
  
  # Console output
  cat("\n========================================================================\n")
  cat("📊 [TRINITY 2/3]: relevant clone fate inferred by isCORAL\n")
  cat("========================================================================\n")
  print(inferred_matrix)
  cat("========================================================================\n\n")
  
  # Color scale
  metaclone_colors <- GetMetacloneColors(n = length(unique_clusters), cluster_names = unique_clusters)
  axis_colors <- metaclone_colors[rev(clustered_order)]
  
  # Dendrogram
  dhc <- ggdendro::dendro_data(hc, type = "rectangle")
  N <- length(unique_clusters)
  p_dendro <- ggplot(ggdendro::segment(dhc)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5, color = "grey40") +
    coord_flip() + scale_x_continuous(limits = c(N + 0.5, 0.5), expand = c(0, 0)) +
    scale_y_reverse(expand = expansion(mult = c(0.1, 0.02))) + theme_void() + theme(plot.margin = margin(0, -5, 0, 0))
  
  # Heatmap
  df_inferred <- reshape2::melt(inferred_matrix)
  colnames(df_inferred) <- c("Metaclone", "CellType", "ExpectedPct")
  df_inferred$Metaclone <- factor(df_inferred$Metaclone, levels = rev(clustered_order))
  df_inferred$CellType <- factor(df_inferred$CellType, levels = target_lineages)
  
  p_main <- ggplot(df_inferred, aes(x = CellType, y = Metaclone, fill = ExpectedPct)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(ExpectedPct > 0.1, ExpectedPct, "")), size = 3.2, color = "black", family = "Arial") +
    scale_fill_gradient(low = "#FDFBFB", high = "#E31A1C", name = "Inferred Fate (%)", limits = c(0, 100)) +
    scale_y_discrete(expand = c(0, 0.5)) +
    theme_minimal(base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
      axis.text.y = element_text(size = 12, color = axis_colors, face = "bold"),
      axis.title = element_text(size = 13), panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 15, face = "plain", margin = margin(b=12)), plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(x = "Inferred Lineage Fate", y = "Metaclone ID", title = "relevant clone fate inferred by isCORAL")
  
  p_combined <- p_dendro + p_main + plot_layout(widths = c(1, 8))
  return(list(table = inferred_matrix, plot = p_combined))
}


#' @title Plot the experimental sibling ground-truth heatmap with a clustering tree
#' @description Calculates the experimental sibling-lineage distribution
#'   associated with each Metaclone.
#' @param seurat_obj A Seurat object containing Metaclone, cell-type, and
#'   experimental clone metadata.
#' @param cluster_col Metadata column containing Metaclone labels.
#' @param cell_type_col Metadata column containing cell-type labels.
#' @param barcode_col Metadata column containing experimental clone labels.
#' @param target_lineages Optional ordered character vector of lineages to display.
#' @param na_barcodes Barcode labels treated as invalid or background.
#' @return A list with the sibling-fate percentage table and combined plot.
#' @export
PlotMetacloneGroundTruthHeatmap <- function(
    seurat_obj, cluster_col = "metaclone", cell_type_col = "cluster_integrated", barcode_col, 
    target_lineages = NULL, na_barcodes = c("NA", "Unknown_NA", "6", NA, "")
) {
  if (missing(barcode_col)) stop("barcode_col is required.")
  
  all_clusters <- as.character(seurat_obj[[cluster_col]][[1]])
  all_cell_types <- as.character(seurat_obj[[cell_type_col]][[1]])
  all_bcs <- as.character(seurat_obj[[barcode_col]][[1]])
  unique_clusters <- as.character(sort(as.numeric(unique(all_clusters))))
  if (is.null(target_lineages)) target_lineages <- sort(unique(all_cell_types))
  
  gt_matrix <- matrix(0, nrow = length(unique_clusters), ncol = length(target_lineages))
  rownames(gt_matrix) <- unique_clusters
  colnames(gt_matrix) <- target_lineages
  
  for (k in unique_clusters) {
    idx <- which(all_clusters == k)
    if(length(idx) == 0) next
    seed_bcs <- all_bcs[idx]
    clean_bcs <- unique(seed_bcs[!(seed_bcs %in% na_barcodes) & !is.na(seed_bcs)])
    
    if(length(clean_bcs) > 0) {
      is_sibling <- (all_bcs %in% clean_bcs)
      is_sibling[idx] <- FALSE # Exclude the seed cluster itself.
      raw_sib_counts <- table(all_cell_types[is_sibling])
      total_sibs <- sum(raw_sib_counts)
      
      if (total_sibs > 0) {
        rel_sib_counts <- (raw_sib_counts / total_sibs) * 100
        for (lin in target_lineages) if (lin %in% names(rel_sib_counts)) gt_matrix[k, lin] <- round(as.numeric(rel_sib_counts[lin]), 1)
      }
    }
  }
  
  # Hierarchical clustering
  hc <- hclust(dist(gt_matrix))
  clustered_order <- rownames(gt_matrix)[hc$order]
  gt_matrix <- gt_matrix[clustered_order, , drop = FALSE]
  
  # Console output
  cat("\n========================================================================\n")
  cat("📊 [TRINITY 3/3]: relevant clone fate groundtruth\n")
  cat("========================================================================\n")
  print(gt_matrix)
  cat("========================================================================\n\n")
  
  # Color scale
  metaclone_colors <- GetMetacloneColors(n = length(unique_clusters), cluster_names = unique_clusters)
  axis_colors <- metaclone_colors[rev(clustered_order)]
  
  # Dendrogram
  dhc <- ggdendro::dendro_data(hc, type = "rectangle")
  N <- length(unique_clusters)
  p_dendro <- ggplot(ggdendro::segment(dhc)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5, color = "grey40") +
    coord_flip() + scale_x_continuous(limits = c(N + 0.5, 0.5), expand = c(0, 0)) +
    scale_y_reverse(expand = expansion(mult = c(0.1, 0.02))) + theme_void() + theme(plot.margin = margin(0, -5, 0, 0))
  
  # Heatmap
  df_gt <- reshape2::melt(gt_matrix)
  colnames(df_gt) <- c("Metaclone", "CellType", "CapturePct")
  df_gt$Metaclone <- factor(df_gt$Metaclone, levels = rev(clustered_order))
  df_gt$CellType <- factor(df_gt$CellType, levels = target_lineages)
  
  p_main <- ggplot(df_gt, aes(x = CellType, y = Metaclone, fill = CapturePct)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(CapturePct > 0.1, CapturePct, "")), size = 3.2, color = "black", family = "Arial") +
    scale_fill_gradient(low = "#F4F9F4", high = "#2E7D32", name = "Sibling Fate (%)", limits = c(0, 100)) +
    scale_y_discrete(expand = c(0, 0.5)) +
    theme_minimal(base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
      axis.text.y = element_text(size = 12, color = axis_colors, face = "bold"),
      axis.title = element_text(size = 13), panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 15, face = "plain", margin = margin(b=12)), plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(x = "Ground Truth Sibling Lineage", y = "Metaclone ID", title = "relevant clone fate groundtruth")
  
  p_combined <- p_dendro + p_main + plot_layout(widths = c(1, 8))
  return(list(table = gt_matrix, plot = p_combined))
}


#' @title Plot the complete inferred-probability facet matrix
#' @description Draws the inferred matrix with seed cells shown in darker cluster colors.
#' @param seurat_obj A Seurat object containing clone-space coordinates,
#'   Metaclone labels, and cached pseudoclone probabilities.
#' @param reduction Reduction used for the inference facets.
#' @param space_title Text identifying the plotted space.
#' @param cluster_col Metadata column containing Metaclone labels.
#' @param custom_colors Optional named vector of cluster colors.
#' @return A patchwork object containing the complete inferred-probability
#'   facet matrix.
#' @export
PlotInferenceSuite <- function(
    seurat_obj, reduction, space_title, 
    cluster_col = "metaclone", 
    custom_colors = NULL
) {
  if (is.null(seurat_obj@misc$isCORAL_embedding$P_norm)) stop("P_norm was not found in seurat_obj@misc$isCORAL_embedding.")
  
  # Extract reduction coordinates and metadata.
  mat <- Seurat::Embeddings(seurat_obj, reduction)
  df_base <- data.frame(
    cell = colnames(seurat_obj), Dim1 = mat[, 1], Dim2 = mat[, 2],
    Cluster = as.character(seurat_obj[[cluster_col]][[1]]),
    stringsAsFactors = FALSE
  )
  
  all_energy_clusters <- seurat_obj[[cluster_col]][[1]]
  P_norm <- seurat_obj@misc$isCORAL_embedding$P_norm
  
  # Configure cluster colors.
  cluster_levels <- as.character(sort(as.numeric(unique(all_energy_clusters))))
  if (is.null(custom_colors)) {
    custom_colors <- GetMetacloneColors(n = length(cluster_levels), cluster_names = cluster_levels)
  }
  
  # ==============================================================================
  # Build and combine the base and darkened seed palettes.
  # ==============================================================================
  darken_color <- function(hex, factor = 0.45) {
    rgb_val <- as.numeric(col2rgb(hex)) * factor
    rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255)
  }
  seed_colors <- sapply(custom_colors, darken_color)
  names(seed_colors) <- paste0(names(custom_colors), "_seed") # For example, "1_seed".
  
  # Combine bright cluster colors with darker seed colors.
  full_palette <- c(custom_colors, seed_colors)
  
  # ==============================================================================
  # Build the plotting data by extracting inferred probabilities per cluster.
  # ==============================================================================
  long_list <- list()
  for (k in sort(unique(all_energy_clusters))) {
    idx <- which(all_energy_clusters == k)
    if (length(idx) == 0) next
    
    prob <- as.numeric(Matrix::colMeans(P_norm[idx, , drop = FALSE]))
    prob[is.na(prob)] <- 0
    if (max(prob) > 0) prob = prob / max(prob)
    
    tmp <- df_base[, c("cell", "Dim1", "Dim2", "Cluster")]
    tmp$TargetCluster <- as.character(k)
    tmp$Prob <- prob
    tmp$IsSeed <- (tmp$Cluster == as.character(k))
    long_list[[length(long_list) + 1]] <- tmp
  }
  
  df_all <- dplyr::bind_rows(long_list)
  df_all$TargetCluster <- factor(df_all$TargetCluster, levels = cluster_levels)
  
  # ==============================================================================
  # Use a dedicated render group to keep facet colors stable.
  # ==============================================================================
  df_all$RenderGroup <- as.character(df_all$TargetCluster)
  df_all$RenderGroup[df_all$IsSeed] <- paste0(df_all$TargetCluster[df_all$IsSeed], "_seed")
  
  # Sort layers so background cells are drawn before high-probability and seed cells.
  df_inferred <- df_all[order(df_all$Prob, df_all$IsSeed), ]
  
  # Base plotting theme.
  base_theme <- theme_classic(base_family = "Arial") +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey95", color = "grey80"),
          strip.text = element_text(size = 11, face = "plain"), axis.text = element_blank(), 
          axis.ticks = element_blank(), axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 18, face = "plain", margin = margin(b=15)))
  
  # ==============================================================================
  # Render inferred probabilities only.
  # ==============================================================================
  p_inferred <- ggplot() +
    geom_point(data = subset(df_inferred, Prob <= 0.02 & !IsSeed), aes(x = Dim1, y = Dim2), color = "#E0E4E7", size = 0.3) +
    # Map the probability cloud to bright RenderGroup colors.
    geom_point(data = subset(df_inferred, Prob > 0.02 & !IsSeed), aes(x = Dim1, y = Dim2, color = RenderGroup, alpha = Prob), size = 0.6) +
    # Map seed cells to darker RenderGroup colors using shape 16.
    geom_point(data = subset(df_inferred, IsSeed), aes(x = Dim1, y = Dim2, color = RenderGroup), shape = 16, size = 1.2, alpha = 1.0) +
    scale_color_manual(values = full_palette) + 
    scale_alpha_continuous(range = c(0, 1.0), guide = "none") + 
    facet_wrap(~ TargetCluster, nrow = 2) + 
    base_theme + 
    labs(title = sprintf("inferred metaclones (%s)", space_title))
  
  return(p_inferred)
}
