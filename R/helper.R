#' @title (v4.0) Calculate Omega-Square for Heritable Genes
#' @description 
#' No mean filtering, consistent with the v2.5 gene set.
#' This version serves as a Seurat pipeline function and invokes the package's internal gene-wise omega-squared engine.
#' Inherits highly safe sparse matrix mechanisms to prevent memory overflow and global chunk-based calculation capabilities.
#' Automatically initializes 'groundtruth', compatible with Seurat V4/V5.
#' On every call, the current normalized data matrix and the requested metadata
#' column are read again from the Seurat object. Legacy cached input stored in
#' \code{@misc$groundtruth$internal_data} is discarded so that changes to
#' \code{group_by}, filtering parameters, the default assay, or expression data
#' are always respected.
#' 
#' @param seurat_obj Seurat object.
#' @param group_by Column name in meta.data used for grouping. MUST be specified (cannot be NULL).
#' @param invalid_barcodes Invalid cell group labels to filter out.
#' @param min_cells_per_group Minimum number of cells required per clone/group. Groups with fewer cells will be filtered out. Default is 2.
#' @param permutation_repeats Number of repeats for the permutation test.
#' @param n_cores Number of cores for parallel computation.
#' @param chunk_size Number of genes per chunk for calculation, controlling memory usage.
#' @param verbose Whether to print progress information.
#' @return The modified Seurat object with observed and permuted omega-squared
#'   results, together with lightweight input provenance, stored in
#'   \code{seurat_obj@misc$groundtruth}.
#' 
#' @export
calculate_heritable_genes <- function(seurat_obj,
                                      group_by = NULL,
                                      invalid_barcodes = c("NA", "Unknown", "unknown", "Unassigned", ""), 
                                      min_cells_per_group = 2,
                                      permutation_repeats = 5,
                                      n_cores = NULL,
                                      chunk_size = 2000,
                                      verbose = TRUE) { 
  
  # --- 1. Initialization and Data Extraction ---
  if (is.null(seurat_obj@misc$groundtruth)) {
    if(verbose) message("Initializing 'groundtruth' in seurat_obj@misc...")
    seurat_obj@misc$groundtruth <- list()
  }
  
  if(verbose) message("Starting Omega-Square calculation using core engine...")

  # Enforce strict check for group_by parameter before resolving the input.
  if (is.null(group_by) || length(group_by) != 1L || is.na(group_by) || !nzchar(group_by)) {
    stop("Error: 'group_by' parameter cannot be NULL. Please specify a valid column name from meta.data.")
  } else if (!(group_by %in% colnames(seurat_obj@meta.data))) {
    stop(paste0("Error: The specified group_by '", group_by, "' was not found in the Seurat object's meta.data."))
  }

  # Never reuse the legacy expression/label cache. It did not encode the
  # parameters or object state that produced it and could therefore become
  # stale when the same object was analysed again with different inputs.
  if (!is.null(seurat_obj@misc$groundtruth$internal_data)) {
    if(verbose) message("   - Removing legacy cached ground-truth input...")
    seurat_obj@misc$groundtruth$internal_data <- NULL
  }

  if(verbose) message("   - Extracting current data directly from Seurat object...")
  assay_used <- SeuratObject::DefaultAssay(seurat_obj)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    sc_expression <- SeuratObject::GetAssayData(
      seurat_obj,
      assay = assay_used,
      layer = "data"
    )
  } else {
    sc_expression <- SeuratObject::GetAssayData(
      seurat_obj,
      assay = assay_used,
      slot = "data"
    )
  }

  expression_cells <- colnames(sc_expression)
  metadata_cells <- rownames(seurat_obj@meta.data)
  if (is.null(expression_cells) || anyNA(expression_cells) ||
      any(!nzchar(expression_cells)) || anyDuplicated(expression_cells)) {
    stop("Error: The current assay data must have unique, non-empty cell names.")
  }
  if (is.null(metadata_cells) || anyNA(metadata_cells) ||
      any(!nzchar(metadata_cells)) || anyDuplicated(metadata_cells)) {
    stop("Error: The Seurat metadata must have unique, non-empty cell names.")
  }

  missing_metadata_cells <- setdiff(expression_cells, metadata_cells)
  if (length(missing_metadata_cells) > 0L) {
    stop(
      paste0(
        "Error: ", length(missing_metadata_cells),
        " assay cell(s) were not found in the Seurat metadata."
      )
    )
  }

  barcode_used <- as.character(
    seurat_obj@meta.data[expression_cells, group_by, drop = TRUE]
  )
  n_cells_input <- ncol(sc_expression)
  
  # --- 2. Remove Invalid Groups & Single-cell Clones ---
  # Mark predefined invalid barcodes
  is_invalid <- is.na(barcode_used) | (barcode_used %in% invalid_barcodes)
  
  # Count cells per clone/group and identify groups smaller than min_cells_per_group
  valid_barcodes_temp <- barcode_used[!is_invalid]
  group_counts <- table(valid_barcodes_temp)
  invalid_clones <- names(group_counts)[group_counts < min_cells_per_group]
  
  if (length(invalid_clones) > 0) {
    if(verbose) message(paste0("   - Identified ", length(invalid_clones), " clone(s) with fewer than ", min_cells_per_group, " cells. Marking for removal..."))
    is_invalid <- is_invalid | (barcode_used %in% invalid_clones)
  }
  
  # Execute filtering operation
  if (any(is_invalid)) {
    if(verbose) message(paste0("   - Removing a total of ", sum(is_invalid), " invalid cells (including small clones)..."))
    valid_idx <- which(!is_invalid)
    sc_expression <- sc_expression[, valid_idx, drop = FALSE]
    barcode_used <- barcode_used[valid_idx]
    
  }
  
  if (is.null(sc_expression) || length(barcode_used) == 0) stop("Error: No valid cells remaining after filtering.")
  
  # --- 3. Core Calculation: OBSERVED Omega-Square ---
  if(verbose) message("   - Step 1/2: Calculating OBSERVED Omega-Square...")
  
  # Directly call the core calculation engine within the package
  omega_square_obs <- CalculateGeneOmegaSquared(
    expression_matrix = sc_expression, 
    groups = barcode_used, 
    chunk_size = chunk_size, 
    verbose = verbose
  )
  
  # Truncate mathematically negative values to 0 for biological relevance
  omega_square_obs_truncated <- pmax(as.numeric(omega_square_obs), 0)
  
  # Assemble the results dataframe
  anova_results <- data.frame(
    name = names(omega_square_obs),
    Omega_square = omega_square_obs_truncated, 
    stringsAsFactors = FALSE
  )
  
  # --- 4. Core Calculation: Permutation Test ---
  if(verbose) message(paste0("   - Step 2/2: Running PERMUTATION (", permutation_repeats, " repeats)..."))
  
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  use_fork <- .Platform$OS.type == "unix" && requireNamespace("doMC", quietly = TRUE)
  
  current_backend <- tryCatch(foreach::getDoParName(), error = function(e) "NONE")
  if (is.null(current_backend) || length(current_backend) == 0) current_backend <- "NONE"
  
  if (use_fork) {
    if(!current_backend %in% c("doMC", "MC")) doMC::registerDoMC(cores = n_cores)
  } else {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    # Safe exit mechanism: Ensure the cluster is properly closed upon errors or user interruption
    on.exit(parallel::stopCluster(cl), add = TRUE) 
  }
  
  opts <- list(preschedule = FALSE)
  
  omega_square_fake_list <- foreach::foreach(
    i = 1:permutation_repeats, 
    .packages = c("Matrix"),
    .export = c("CalculateGeneOmegaSquared"), 
    .options.multicore = opts
  ) %dopar% {
    
    # Only shuffle labels
    groups_fake <- sample(barcode_used)
    
    fake_res <- CalculateGeneOmegaSquared(
      expression_matrix = sc_expression, 
      groups = groups_fake, 
      chunk_size = chunk_size, 
      verbose = FALSE 
    )
    
    # Truncate background noise output to ensure valid Threshold calculation
    return(pmax(as.numeric(fake_res), 0))
  }
  
  omega_square_fake <- unlist(omega_square_fake_list)
  threshold <- 1.5 * stats::IQR(omega_square_fake, na.rm = TRUE) + stats::quantile(omega_square_fake, 0.75, na.rm = TRUE)
  if(verbose) message(paste0("   - Statistical threshold calculated: ", round(threshold, 5)))
  
  # --- 5. Write Results ---
  seurat_obj@misc$groundtruth$heritable_genes_df <- anova_results
  seurat_obj@misc$groundtruth$heritable_genes_threshold <- threshold 
  seurat_obj@misc$groundtruth$omega_squared_null_distribution <- omega_square_fake
  seurat_obj@misc$groundtruth$input_provenance <- list(
    group_by = group_by,
    invalid_barcodes = invalid_barcodes,
    min_cells_per_group = min_cells_per_group,
    assay = assay_used,
    layer = "data",
    n_cells_input = n_cells_input,
    n_cells_used = ncol(sc_expression),
    n_cells_removed = n_cells_input - ncol(sc_expression),
    n_groups_used = length(unique(barcode_used)),
    group_sizes_used = table(barcode_used),
    groups_removed_for_size = invalid_clones,
    permutation_repeats = permutation_repeats,
    n_cores = n_cores,
    chunk_size = chunk_size,
    seurat_object_version = as.character(utils::packageVersion("SeuratObject"))
  )
  
  if(verbose) message("Omega-Square pipeline complete.")
  return(seurat_obj)
}
