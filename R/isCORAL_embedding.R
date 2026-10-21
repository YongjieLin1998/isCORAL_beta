#' @title Unified Integrated isCORAL Dual-Layer Manifold Embedding Engine (CORALembedding)
#' @description 
#' Directly fits the latent space of pseudoclones using Seurat-normalized transcriptomic features,
#' then applies full-weight graph diffusion at the output space.
#' Permanently mounts "CORAL_cell_embed" (pure gene cell space) and "CORAL_clone_embed" (fully diffused clone space),
#' while automatically caching the generated normalized graph probability, MDS scaffold, and initial projections 
#' into @misc$isCORAL_embedding for persistent reuse.
#' 
#' @param seurat_obj A preprocessed (Normalize/Scale) Seurat object compatible with V4 and V5.
#' @param binsize_key Character key for locating the graph transition matrix inside @misc. Default is "binsize_200".
#' @param n_dims Target dimensionality of the latent manifold space (typically 20 or 30).
#' @param top_n_genes Number of core driving genes filtered by inferred
#'   Omega-Square ranking. At least two ranked genes must also be present in
#'   the normalized RNA data.
#' @param max_cells_lasso Maximum number of cells subsampled for training the regression mapper. 
#'        Default is `NULL` (subsampling disabled, uses all cells). Set to an integer (e.g., 3000) to enable subsampling for very large datasets to prevent memory exhaustion. The resulting training set must contain at least five cells for five-fold cross-validation.
#' @param alpha_elasticnet Elastic net penalty coefficient. Default is 0.05.
#' @param workers_glmnet Number of parallel threads used during cross-validation training, hard-capped at 5 internally.
#' @param random_seed Optional integer seed for the MDS, training-cell
#'   subsampling, and cross-validation steps. The default is \code{NULL}, which
#'   does not reset the caller's random-number generator and preserves the
#'   stochastic behavior of version 6.1.4.
#' @param verbose Whether to print routine progress messages. Defaults to
#'   \code{FALSE}; warnings and errors are not suppressed.
#' 
#' @return A Seurat object containing two standard dimensional reductions and
#'   a backward-compatible persistent cache in
#'   \code{@misc$isCORAL_embedding}. Beta 1.0.0 retains all legacy cache
#'   fields and additionally stores requested and effective parameters,
#'   the actual normalized-data access path and selection strategy, declared
#'   dependency defaults, random-number states, training cells, data
#'   provenance, and software versions.
#' @importFrom stats predict
#' @export
CORALembedding <- function(
    seurat_obj, 
    binsize_key = "binsize_200", 
    n_dims = 30, 
    top_n_genes = 3000, 
    max_cells_lasso = NULL,          
    alpha_elasticnet = 0.05, 
    workers_glmnet = 5,
    random_seed = NULL,
    verbose = FALSE
) {
  embedding_call <- match.call()

  capture_rng_state <- function() {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      return(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
    NULL
  }

  package_version_or_na <- function(package_name) {
    tryCatch(
      as.character(utils::packageVersion(package_name)),
      error = function(e) NA_character_
    )
  }

  if (!is.character(binsize_key) || length(binsize_key) != 1 ||
      is.na(binsize_key) || !nzchar(binsize_key)) {
    stop("'binsize_key' must be one non-empty character value.")
  }
  if (!is.numeric(n_dims) || length(n_dims) != 1 || is.na(n_dims) ||
      !is.finite(n_dims) || n_dims < 2 || n_dims != floor(n_dims)) {
    stop("'n_dims' must be one integer greater than or equal to 2.")
  }
  if (!is.numeric(top_n_genes) || length(top_n_genes) != 1 ||
      is.na(top_n_genes) || !is.finite(top_n_genes) ||
      top_n_genes < 1 || top_n_genes != floor(top_n_genes)) {
    stop("'top_n_genes' must be one positive integer.")
  }
  if (!is.numeric(alpha_elasticnet) || length(alpha_elasticnet) != 1 ||
      is.na(alpha_elasticnet) || !is.finite(alpha_elasticnet) ||
      alpha_elasticnet < 0 || alpha_elasticnet > 1) {
    stop("'alpha_elasticnet' must be one number in the closed interval [0, 1].")
  }
  if (!is.numeric(workers_glmnet) || length(workers_glmnet) != 1 ||
      is.na(workers_glmnet) || !is.finite(workers_glmnet) ||
      workers_glmnet < 1 || workers_glmnet != floor(workers_glmnet)) {
    stop("'workers_glmnet' must be one positive integer.")
  }
  if (!is.null(max_cells_lasso) &&
      (!is.numeric(max_cells_lasso) || length(max_cells_lasso) != 1 ||
       is.na(max_cells_lasso) ||
       (!is.infinite(max_cells_lasso) &&
        (max_cells_lasso < 1 || max_cells_lasso != floor(max_cells_lasso))) ||
       (is.infinite(max_cells_lasso) && max_cells_lasso < 0))) {
    stop("'max_cells_lasso' must be NULL, Inf, or one positive integer.")
  }
  if (!is.null(random_seed) &&
      (!is.numeric(random_seed) || length(random_seed) != 1 ||
       is.na(random_seed) || !is.finite(random_seed) ||
       random_seed < 0 || random_seed > .Machine$integer.max ||
       random_seed != floor(random_seed))) {
    stop("'random_seed' must be NULL or one non-negative integer.")
  }

  n_dims <- as.integer(n_dims)
  top_n_genes <- as.integer(top_n_genes)
  workers_glmnet <- as.integer(workers_glmnet)
  if (!is.null(max_cells_lasso) && is.finite(max_cells_lasso)) {
    max_cells_lasso <- as.integer(max_cells_lasso)
  }
  if (!is.null(random_seed)) random_seed <- as.integer(random_seed)

  rng_kind_at_entry <- RNGkind()
  rng_state_at_entry <- capture_rng_state()
  if (!is.null(random_seed)) set.seed(random_seed)
  rng_state_after_seed <- capture_rng_state()

  old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  Sys.setenv(OMP_NUM_THREADS = "1")
  on.exit(.isCORAL_restore_omp_threads(old_omp), add = TRUE)
  
  if (verbose) message("🚀 [isCORAL] Initiating unified manifold embedding pipeline...")
  
  coral_misc_key <- ifelse(
    "CORAL_scan_results" %in% names(seurat_obj@misc),
    "CORAL_scan_results",
    "coral_scan_results"
  )
  coral_res <- seurat_obj@misc[[coral_misc_key]][[binsize_key]]
  if (is.null(coral_res)) stop(sprintf("🚨 Fatal Error: Unable to locate network scan results for key '%s' inside @misc!", binsize_key))
  
  prob_matrix <- coral_res$S_stab
  ss_stats <- coral_res$ss_filter_average
  if (is.null(prob_matrix) || is.null(ss_stats) ||
      !"omega_square_infer" %in% colnames(ss_stats)) {
    stop(sprintf(
      "🚨 Fatal Error: The scan result '%s' does not contain S_stab and a valid ss_filter_average table.",
      binsize_key
    ))
  }
  
  read_normalized_layer <- function() {
    list(
      matrix = Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data"),
      api = "layer"
    )
  }
  read_normalized_slot <- function() {
    list(
      matrix = Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "data"),
      api = "slot"
    )
  }
  seurat_object_version <- tryCatch(
    utils::packageVersion("SeuratObject"),
    error = function(e) NULL
  )
  layer_first <- !is.null(seurat_object_version) &&
    seurat_object_version >= "5.0.0"
  expression_access <- if (layer_first) {
    tryCatch(
      read_normalized_layer(),
      error = function(e) read_normalized_slot()
    )
  } else {
    tryCatch(
      read_normalized_slot(),
      error = function(e) read_normalized_layer()
    )
  }
  E_mat <- expression_access$matrix
  expression_access_api <- expression_access$api
  expression_access_strategy <- if (layer_first) {
    "SeuratObject >= 5: layer first, slot fallback"
  } else {
    "SeuratObject < 5 or unavailable: slot first, layer fallback"
  }
  expression_matrix_class <- class(E_mat)
  expression_matrix_dimensions <- dim(E_mat)
  rm(expression_access)
  
  seurat_cells <- colnames(E_mat)
  N_cells_raw <- length(seurat_cells)
  N_genes_raw <- nrow(E_mat)
  input_object_descriptor <- list(
    storage = "descriptor_only",
    class = class(seurat_obj),
    project_name = tryCatch(seurat_obj@project.name, error = function(e) NA_character_),
    n_cells = ncol(seurat_obj),
    n_features = nrow(seurat_obj),
    default_assay = Seurat::DefaultAssay(seurat_obj)
  )
  
  if (length(dim(prob_matrix)) != 2L || nrow(prob_matrix) != ncol(prob_matrix)) {
    stop("🚨 The probabilistic matrix must be a square two-dimensional matrix.")
  }

  probability_dimnames_assigned <- FALSE
  probability_dimname_action <- "none"
  probability_rownames <- rownames(prob_matrix)
  probability_colnames <- colnames(prob_matrix)

  if (is.null(probability_rownames) && is.null(probability_colnames)) {
    if (nrow(prob_matrix) != N_cells_raw) {
      stop("🚨 Dimension mismatch between expression matrix and probabilistic matrix! Please align cell barcodes beforehand.")
    }
    rownames(prob_matrix) <- colnames(prob_matrix) <- seurat_cells
    probability_dimnames_assigned <- TRUE
    probability_dimname_action <- "assigned_from_seurat_cells"
  } else if (is.null(probability_rownames)) {
    rownames(prob_matrix) <- probability_colnames
    probability_dimnames_assigned <- TRUE
    probability_dimname_action <- "copied_from_columns"
  } else if (is.null(probability_colnames)) {
    colnames(prob_matrix) <- probability_rownames
    probability_dimnames_assigned <- TRUE
    probability_dimname_action <- "copied_from_rows"
  }

  probability_rownames <- rownames(prob_matrix)
  probability_colnames <- colnames(prob_matrix)
  invalid_probability_names <- function(cell_names) {
    is.null(cell_names) || length(cell_names) != nrow(prob_matrix) ||
      anyNA(cell_names) || any(!nzchar(cell_names)) || anyDuplicated(cell_names)
  }
  if (invalid_probability_names(probability_rownames) ||
      invalid_probability_names(probability_colnames)) {
    stop("🚨 The probabilistic matrix must have complete, unique cell barcodes on both axes.")
  }
  if (!setequal(probability_rownames, probability_colnames)) {
    stop("🚨 The row and column cell-barcode sets of the probabilistic matrix must be identical.")
  }
  
  valid_cells <- intersect(
    seurat_cells,
    intersect(probability_rownames, probability_colnames)
  )
  N_cells <- length(valid_cells)
  if (N_cells == 0) stop("🚨 Fatal Error: No valid overlapping cells found between expression matrix and prob_matrix.")
  
  ss_stats_sorted <- ss_stats[order(ss_stats$omega_square_infer, decreasing = TRUE), , drop = FALSE]
  top_genes <- intersect(rownames(ss_stats_sorted)[1:min(top_n_genes, nrow(ss_stats_sorted))], rownames(E_mat))
  if (length(top_genes) < 2L) {
    stop("🚨 Fatal Error: At least two ranked memory genes must be present in the RNA data layer for glmnet fitting.")
  }
  
  E_mat <- E_mat[top_genes, valid_cells, drop = FALSE]
  
  P <- prob_matrix[valid_cells, valid_cells, drop = FALSE]
  if (!inherits(P, "CsparseMatrix")) {
    P <- as(P, "CsparseMatrix")
  }
  
  row_sums <- Matrix::rowSums(P)
  row_sums[row_sums == 0] <- 1  
  P_norm <- Matrix::Diagonal(x = 1 / row_sums) %*% P
  
  v_norms <- Matrix::colSums(E_mat^2)
  chunk_size <- ifelse(N_cells < 15000, 5000, 2000) 
  starts <- seq(1, N_cells, by = chunk_size)
  E_cross <- matrix(0, nrow = N_cells, ncol = N_cells)
  
  for (s in starts) {
    idx <- s:min(s + chunk_size - 1, N_cells)
    cross_chunk <- as.matrix(Matrix::crossprod(E_mat, E_mat[, idx, drop = FALSE])) 
    D_sq_chunk <- -2 * cross_chunk + v_norms + rep(v_norms[idx], each = N_cells)
    D_sq_chunk[D_sq_chunk < 0] <- 0 
    E_cross <- E_cross + as.matrix(P_norm %*% sqrt(D_sq_chunk)) %*% as.matrix(t(P_norm[, idx, drop = FALSE]))
  }
  
  E_in <- diag(E_cross)
  Dist_E_MMD2 <- 2 * E_cross
  Dist_E_MMD2 <- sweep(Dist_E_MMD2, 1, E_in, "-")
  Dist_E_MMD2 <- sweep(Dist_E_MMD2, 2, E_in, "-")
  Dist_E_MMD2[Dist_E_MMD2 < 0] <- 0 
  rm(E_cross)
  
  row_mean <- rowMeans(Dist_E_MMD2)
  mean_all <- mean(row_mean)
  B_mat <- -0.5 * (sweep(sweep(Dist_E_MMD2, 1, row_mean, "-"), 2, row_mean, "-") + mean_all)
  
  safe_n_dims <- as.integer(min(n_dims, N_cells - 5))
  if (safe_n_dims < 2) stop("🚨 Not enough cells to perform dimensional reduction. N_cells too small.")
  
  rng_state_before_irlba <- capture_rng_state()
  svd_res <- irlba::irlba(B_mat, nv = safe_n_dims, fastpath = TRUE, tol = 1e-4)
  rng_state_after_irlba <- capture_rng_state()
  
  original_MDS_coords <- sweep(svd_res$u, 2, sqrt(pmax(svd_res$d[1:safe_n_dims], 0)), "*")
  rownames(original_MDS_coords) <- valid_cells
  
  rm(B_mat, Dist_E_MMD2, svd_res); gc(verbose = FALSE)
  
  X_raw <- t(as.matrix(E_mat))
  train_cell_norms <- sqrt(rowSums(X_raw^2) + 1e-8)
  X_l2 <- sweep(X_raw, 1, train_cell_norms, "/")
  ref_gene_means <- matrixStats::colMeans2(X_l2)
  ref_gene_sds <- matrixStats::colSds(X_l2)
  ref_gene_sds[ref_gene_sds == 0] <- 1e-8
  X_scaled <- as.matrix(scale(X_l2, center = ref_gene_means, scale = ref_gene_sds))
  X_scaled[is.na(X_scaled)] <- 0
  
  do_subsample <- TRUE
  if (is.null(max_cells_lasso) || is.infinite(max_cells_lasso)) {
    do_subsample <- FALSE
  } else if (nrow(X_scaled) <= max_cells_lasso) {
    do_subsample <- FALSE
  }
  
  rng_state_before_subsampling <- capture_rng_state()
  if (do_subsample) {
    train_idx <- base::sample(1:nrow(X_scaled), max_cells_lasso)
    if (verbose) message(sprintf("📉 [isCORAL] Subsampling enabled: Randomly selected %d cells for ElasticNet training.", max_cells_lasso))
  } else {
    train_idx <- 1:nrow(X_scaled)
    if (verbose) message(sprintf("📈 [isCORAL] Subsampling bypassed: Using all %d cells for ElasticNet training.", length(train_idx)))
  }
  rng_state_after_subsampling <- capture_rng_state()
  training_cells <- valid_cells[train_idx]

  glmnet_nfolds <- 5L
  if (length(train_idx) < glmnet_nfolds) {
    stop(sprintf(
      "Five-fold Elastic Net cross-validation requires at least %d training cells; only %d were selected.",
      glmnet_nfolds, length(train_idx)
    ))
  }

  glmnet_control_before_fit <- glmnet::glmnet.control()
  declared_defaults <- list(
    irlba = as.list(formals(irlba::irlba)),
    cv_glmnet = as.list(formals(glmnet::cv.glmnet)),
    glmnet = as.list(formals(glmnet::glmnet))
  )

  describe_foreach_backend <- function() {
    list(
      registered = foreach::getDoParRegistered(),
      name = tryCatch(foreach::getDoParName(), error = function(e) NULL),
      version = tryCatch(foreach::getDoParVersion(), error = function(e) NULL),
      workers = tryCatch(foreach::getDoParWorkers(), error = function(e) NA_integer_)
    )
  }
  foreach_backend_before <- describe_foreach_backend()
  foreach_globals <- get(
    ".foreachGlobals",
    envir = asNamespace("foreach"),
    inherits = FALSE
  )
  previous_foreach_backend <- NULL
  if (isTRUE(foreach_backend_before$registered)) {
    previous_foreach_backend <- list(
      fun = get("fun", envir = foreach_globals, inherits = FALSE),
      data = get("data", envir = foreach_globals, inherits = FALSE),
      info = get("info", envir = foreach_globals, inherits = FALSE)
    )
  }

  detected_cores <- suppressWarnings(parallel::detectCores())
  available_workers <- if (
    length(detected_cores) != 1 || is.na(detected_cores) || detected_cores < 2
  ) 1L else as.integer(detected_cores - 1L)
  actual_workers <- max(1L, min(workers_glmnet, 5L, available_workers))
  cl <- parallel::makeCluster(actual_workers)
  doParallel::registerDoParallel(cl)
  foreach_backend_during <- describe_foreach_backend()
  
  rng_state_before_cv <- capture_rng_state()
  cv_foldid <- base::sample(rep(seq_len(glmnet_nfolds), length.out = length(train_idx)))
  names(cv_foldid) <- training_cells

  cv_fit <- tryCatch({
    glmnet::cv.glmnet(x = X_scaled[train_idx, ], y = original_MDS_coords[train_idx, ], 
                      family = "mgaussian", alpha = alpha_elasticnet,
                      nfolds = glmnet_nfolds, foldid = unname(cv_foldid),
                      parallel = TRUE, standardize = TRUE, intercept = TRUE,
                      standardize.response = FALSE)
  }, finally = {
    try(parallel::stopCluster(cl), silent = TRUE)
    if (is.null(previous_foreach_backend)) {
      tryCatch({
        temporary_backend_bindings <- intersect(
          c("fun", "data", "info"),
          ls(envir = foreach_globals, all.names = TRUE)
        )
        if (length(temporary_backend_bindings) > 0L) {
          rm(list = temporary_backend_bindings, envir = foreach_globals)
        }
      }, error = function(e) foreach::registerDoSEQ())
    } else {
      tryCatch(
        foreach::setDoPar(
          previous_foreach_backend$fun,
          previous_foreach_backend$data,
          previous_foreach_backend$info
        ),
        error = function(e) foreach::registerDoSEQ()
      )
    }
  })
  rng_state_after_cv <- capture_rng_state()
  foreach_backend_after <- describe_foreach_backend()
  glmnet_nfolds_effective <- length(unique(cv_foldid))

  Y_ref_raw <- stats::predict(cv_fit, newx = X_scaled, s = "lambda.min")[, , 1]
  target_dimension_sds <- matrixStats::colSds(original_MDS_coords)
  raw_projection_sds <- matrixStats::colSds(Y_ref_raw)
  variance_scale_factors <- target_dimension_sds / pmax(raw_projection_sds, 1e-8)
  Y_ref_raw <- sweep(Y_ref_raw, 2, variance_scale_factors, "*")
  rownames(Y_ref_raw) <- valid_cells
  colnames(Y_ref_raw) <- paste0("CORALcell_", 1:safe_n_dims)
  
  seurat_obj[["CORAL_cell_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = Y_ref_raw, 
    key = "CORALcell_", 
    assay = Seurat::DefaultAssay(seurat_obj)
  )
  
  Y_smoothed <- as.matrix(P_norm %*% Y_ref_raw)
  rownames(Y_smoothed) <- valid_cells
  colnames(Y_smoothed) <- paste0("CORALclone_", 1:safe_n_dims)
  
  seurat_obj[["CORAL_clone_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = Y_smoothed, 
    key = "CORALclone_", 
    assay = Seurat::DefaultAssay(seurat_obj)
  )
  
  rng_state_at_cache_write <- capture_rng_state()
  run_timestamp_utc <- format(
    as.POSIXct(Sys.time(), tz = "UTC"),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )
  software_versions <- list(
    R = R.version.string,
    platform = R.version$platform,
    Seurat = package_version_or_na("Seurat"),
    SeuratObject = package_version_or_na("SeuratObject"),
    Matrix = package_version_or_na("Matrix"),
    matrixStats = package_version_or_na("matrixStats"),
    irlba = package_version_or_na("irlba"),
    glmnet = package_version_or_na("glmnet"),
    doParallel = package_version_or_na("doParallel"),
    foreach = package_version_or_na("foreach")
  )
  implementation_defaults <- list(
    declared_defaults = declared_defaults,
    glmnet_control_before_fit = glmnet_control_before_fit
  )

  seurat_obj@misc$isCORAL_embedding <- list(
    P_norm = P_norm,                           
    original_MDS_coords = original_MDS_coords,   
    Y_ref_raw = Y_ref_raw,                     
    valid_cells = valid_cells,                 
    top_genes = top_genes,                     
    source_binsize = binsize_key,              
    lasso_model = cv_fit,
    schema_version = "1.0",
    iscoral_version = "1.0.0",
    release_status = "beta",
    method = "pCORAL dual-layer manifold embedding",
    call = embedding_call,
    call_text = paste(deparse(embedding_call), collapse = " "),
    run_parameters = list(
      seurat_obj = input_object_descriptor,
      binsize_key = binsize_key,
      n_dims = n_dims,
      top_n_genes = top_n_genes,
      max_cells_lasso = max_cells_lasso,
      alpha_elasticnet = alpha_elasticnet,
      workers_glmnet = workers_glmnet,
      random_seed = random_seed,
      verbose = verbose
    ),
    effective_parameters = list(
      scan_misc_key = coral_misc_key,
      assay = "RNA",
      normalized_data_name = "data",
      assay_access_api = expression_access_api,
      assay_access_strategy = expression_access_strategy,
      default_assay_for_reductions = Seurat::DefaultAssay(seurat_obj),
      n_input_cells = N_cells_raw,
      n_input_genes = N_genes_raw,
      n_valid_cells = N_cells,
      excluded_cells = setdiff(seurat_cells, valid_cells),
      n_selected_genes = length(top_genes),
      requested_dimensions = n_dims,
      effective_dimensions = safe_n_dims,
      distance_chunk_size = chunk_size,
      subsampling_applied = do_subsample,
      n_training_cells = length(train_idx),
      detected_cores = detected_cores,
      available_workers = available_workers,
      actual_workers = actual_workers,
      irlba = list(nv = safe_n_dims, fastpath = TRUE, tol = 1e-4),
      glmnet = list(
        family = "mgaussian",
        alpha = alpha_elasticnet,
        nfolds = glmnet_nfolds,
        nfolds_requested = glmnet_nfolds,
        nfolds_effective = glmnet_nfolds_effective,
        foldid_supplied = TRUE,
        n_training_observations = length(train_idx),
        n_predictors = ncol(X_scaled),
        n_response_dimensions = ncol(original_MDS_coords),
        type_measure_requested = "default",
        type_measure_effective = cv_fit$name,
        alignment = "lambda",
        grouped_requested = TRUE,
        grouped_effective = length(train_idx) / glmnet_nfolds >= 3,
        weights_requested = NULL,
        weights_effective = "equal weights",
        offset = NULL,
        lambda_requested = NULL,
        lambda_sequence = unname(cv_fit$lambda),
        lambda_path_length = length(cv_fit$lambda),
        lambda_path_ratio = min(cv_fit$lambda) / max(cv_fit$lambda),
        parallel = TRUE,
        keep = FALSE,
        relax = FALSE,
        trace_it = 0,
        standardize = TRUE,
        intercept = TRUE,
        standardize_response = FALSE,
        lambda_selection = "lambda.min",
        lambda_min = unname(cv_fit$lambda.min),
        lambda_1se = unname(cv_fit$lambda.1se),
        cv_fit_call = cv_fit$call,
        glmnet_fit_call = cv_fit$glmnet.fit$call,
        arguments_using_glmnet_declared_defaults = c(
          "nlambda", "lambda.min.ratio", "penalty.factor", "exclude",
          "lower.limits", "upper.limits"
        )
      )
    ),
    reproducibility = list(
      random_seed = random_seed,
      seed_set_by_function = !is.null(random_seed),
      rng_kind_at_entry = rng_kind_at_entry,
      rng_kind_at_cache_write = RNGkind(),
      rng_state = list(
        at_entry = rng_state_at_entry,
        after_seed = rng_state_after_seed,
        before_irlba = rng_state_before_irlba,
        after_irlba = rng_state_after_irlba,
        before_subsampling = rng_state_before_subsampling,
        after_subsampling = rng_state_after_subsampling,
        before_cv = rng_state_before_cv,
        after_cv = rng_state_after_cv,
        at_cache_write = rng_state_at_cache_write
      ),
      training_indices = train_idx,
      training_cells = training_cells,
      cv_foldid = cv_foldid
    ),
    preprocessing = list(
      cell_l2_normalization = TRUE,
      l2_epsilon = 1e-8,
      gene_z_score_scaling = TRUE,
      reference_gene_means = ref_gene_means,
      reference_gene_sds = ref_gene_sds,
      target_dimension_sds = target_dimension_sds,
      raw_projection_sds = raw_projection_sds,
      variance_scale_factors = variance_scale_factors
    ),
    provenance = list(
      timestamp_utc = run_timestamp_utc,
      scan_misc_key = coral_misc_key,
      source_binsize = binsize_key,
      probability_matrix_class = class(prob_matrix),
      probability_matrix_dimensions = dim(prob_matrix),
      probability_dimnames_assigned = probability_dimnames_assigned,
      probability_dimname_action = probability_dimname_action,
      expression_matrix_class = expression_matrix_class,
      expression_matrix_dimensions = expression_matrix_dimensions,
      software_versions = software_versions,
      external_software = as.list(extSoftVersion()),
      session_info = utils::sessionInfo(),
      implementation_defaults = implementation_defaults,
      foreach_backend = list(
        before = foreach_backend_before,
        during = foreach_backend_during,
        after = foreach_backend_after
      ),
      omp_num_threads_before = old_omp,
      omp_num_threads_during = "1",
      omp_num_threads_after = .isCORAL_normalize_omp_threads(old_omp)
    ),
    outputs = list(
      cell_reduction = "CORAL_cell_embed",
      cell_reduction_key = "CORALcell_",
      clone_reduction = "CORAL_clone_embed",
      clone_reduction_key = "CORALclone_",
      assay = Seurat::DefaultAssay(seurat_obj),
      n_cells = N_cells,
      dimensions = safe_n_dims
    )
  )
  
  return(seurat_obj)
}

#' @title Dynamic manifold spectral smoothing
#' @description Recomputes the clone-space embedding from the persistent
#'   isCORAL embedding cache at a selected graph-smoothing strength.
#' @param seurat_obj A Seurat object produced by \code{CORALembedding()}.
#' @param lambda Smoothing coefficient in the closed interval from 0 to 1.
#' @return The modified Seurat object with a new reduction whose name begins
#'   with \code{CORAL_smooth_lambda_}.
#' @export
SmoothEmbedding <- function(seurat_obj, lambda = 0.5) {
  
  if (is.null(seurat_obj@misc$isCORAL_embedding)) {
    stop("🚨 Error: 'isCORAL_embedding' persistent cache not detected in the object! Please run CORALembedding() first.")
  }
  
  if (lambda < 0 || lambda > 1) {
    stop("🚨 Parameter Out of Bounds: lambda smoothing coefficient must be strictly restricted within the [0, 1] interval!")
  }
  
  cache       <- seurat_obj@misc$isCORAL_embedding
  P_norm      <- cache$P_norm
  Y_ref_raw   <- cache$Y_ref_raw
  valid_cells <- cache$valid_cells
  n_dims      <- ncol(Y_ref_raw)
  
  message(sprintf("⚙️ [isCORAL] Invoking persistent basis, dynamically compiling manifold spectrum for lambda = %s...", lambda))
  
  Y_smooth_max <- as.matrix(P_norm %*% Y_ref_raw)
  Y_final      <- (1 - lambda) * Y_ref_raw + lambda * Y_smooth_max
  
  lambda_tag     <- gsub("\\.", "", as.character(lambda))
  reduction_name <- paste0("CORAL_smooth_lambda_", lambda_tag)
  reduction_key  <- paste0("CORALsmooth", lambda_tag, "_")
  
  rownames(Y_final) <- valid_cells
  colnames(Y_final) <- paste0(reduction_key, 1:n_dims)
  
  seurat_obj[[reduction_name]] <- Seurat::CreateDimReducObject(
    embeddings = Y_final,
    key = reduction_key,
    assay = Seurat::DefaultAssay(seurat_obj)
  )
  return(seurat_obj)
}
