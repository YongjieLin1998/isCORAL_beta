# isCORAL beta 1.0.0

## Release identity

- Beta 1.0.0 packages the source-audited isCORAL 6.1.5 implementation under the
  public package version `1.0.0`.
- The canonical core entry point is now `isCORAL()`. Development-stage entry
  points containing internal version numbers were removed and are not kept as
  aliases.
- Core scan defaults, input modes, Seurat v4/v5 access, pCORAL embedding,
  Metaclone inference, plotting functions, and compatibility wrappers for
  renamed non-core functions are inherited from the 6.1.5 implementation.

## Documentation

- Added the exported `PlotMemoryGenes()` API for resolving scan-result tables,
  filtering genes by mean expression, ranking them by inferred omega-squared,
  and drawing a compact bar chart. `return_data = TRUE` also exposes the full
  filtered table without modifying the stored scan result.
- Added the exported `CompareMemoryScores()` API for aligning isCORAL-inferred
  memory scores with experimental barcode-derived omega-squared values and
  returning a filtered gene table, ROC/PR curves, AUROC, AUPRC, AUPRC
  enrichment, maximum F1, correlations, and RMSE without modifying the Seurat
  object.
- Added `Assign_deterministic_pseudoclone()` as the canonical deterministic
  hierarchical assignment API. It preserves the established `1 - S_stab`,
  `hclust()`, and `cutree()` path, and writes assignments to
  `pseudoclone_k<k>_bs<bin_size>` metadata columns.
- Added `simulate_clonal_expression()` by promoting the isCORAL-only portion of
  the clone-expression benchmark visualization into the package. It accepts
  explicit genes, assay/layer and reduction settings, draws observed expression
  on the original UMAP together with pseudoclone-wise violin distributions,
  and can return the underlying cell-level and summary tables.
- Added a complete UTF-8 Chinese-English version guide covering versions
  6.1.0 through 6.1.5, including scientific defaults, input modes, public API
  availability, pseudoclone and Metaclone terminology, migration profiles,
  reproducibility boundaries, dependencies, and known visualization limits.

## Compatibility

- `leidenbase` remains an optional `Suggests` dependency and is the supported
  Leiden backend for `FindMetaclone()`. The function never installs software
  at runtime; if `algorithm = 4` is requested without an available
  `leidenbase` namespace, or with Seurat < 5.2.0, it warns, uses Louvain
  (`algorithm = 1`), and records the effective fallback. The ambiguous
  Python-backed `leiden` branch was removed.
- `Assign_deterministic_pseudoclone()` is now the canonical name for the simple
  hierarchical clustering of `1 - S_stab`. `SimulateClone()`,
  `Simulate_clonal_expression_simple()`, and `SimulateMetaClone()` are retained
  only as deprecated forwarding aliases. The clustering path, defaults, and
  barcode alignment are unchanged. Output metadata now uses the
  `pseudoclone_k<k>_bs<bin_size>` prefix instead of the legacy
  `metaclone_k<k>_bs<bin_size>` prefix, separating deterministic pseudoclones
  from embedding-based Metaclones.
- `PlotPseudoclone()` now treats `seed_cell` and `random_seed` as the canonical
  parameter names. The older `seed_cells` and `seed` aliases remain accepted
  for migration but emit deprecation warnings.
- `FindMetaclone()` remains the canonical embedding-based Metaclone
  community-detection function, with the same arguments, defaults, graph
  names, metadata behavior, and clustering path as version 6.1.4.
- `FindMetaclones()` is restored only as a deprecated forwarding wrapper for
  scripts written for versions 6.1.0 and 6.1.1. It does not contain a second
  implementation.
- When the historical default `dims = 1:30` is used with a lower-dimensional
  CORAL reduction, `FindMetaclone()` warns and uses the available dimensions.
  Other explicit out-of-range dimension requests remain errors.
- The broad `RcppML` namespace import was removed. `RcppML::nmf()` remains an
  explicit dependency and call, while the previous `predict` and `head`
  masking warnings are avoided.

## Reproducibility and provenance

- `calculate_heritable_genes()` no longer reuses or stores the legacy
  `@misc$groundtruth$internal_data` expression/label cache. Every invocation
  resolves the current default-assay `data` matrix and requested metadata
  labels again, aligns them by cell barcode, and records lightweight input
  provenance. Re-running the function on one object with a different
  `group_by` or filtering configuration therefore cannot silently reuse stale
  inputs.
- `CORALembedding()` retains its seven legacy
  `@misc$isCORAL_embedding` fields without renaming or moving them.
- The cache now records requested and effective parameters, version-guided
  SeuratObject 5.x layer or 4.x slot access, including fallback strategy and
  the actual path used, input and training cell counts,
  selected genes, effective dimensions, actual workers, IRLBA settings,
  Elastic Net settings, selected lambdas, training cells, cross-validation
  folds, random-number states, scaling statistics, reduction names, and
  software versions.
- A final optional `random_seed = NULL` argument was added.
  `NULL` preserves the inherited entry-RNG semantics by not resetting the
  caller's generator; supplying an integer initializes the embedding workflow
  reproducibly.
- `max_cells_lasso = NULL` is recorded explicitly as no subsampling, together
  with `subsampling_applied = FALSE` and the actual number of training cells.
- Five-fold fitting now fails early with a clear error if fewer than five
  training cells are selected, and glmnet fitting requires at least two
  matching memory genes.
- Both axes of the probability matrix are validated by cell barcode. A single
  missing dimname axis is restored from the named axis, while inconsistent
  row and column cell sets are rejected.
- The temporary doParallel backend used by `CORALembedding()` is closed and
  the previous registered or unregistered foreach state is restored; a safe
  sequential backend is used only if restoration fails.
- The cache additionally records the fitted glmnet calls, declared dependency
  defaults, the pre-fit `glmnet.control()` state, foreach backend descriptors,
  external-software information, and `sessionInfo()`.

## Metaclone records

- `FindMetaclone()` continues to write its labels to the selected Seurat
  metadata column and set that column as the active identity.
- Each run now records requested and effective clustering parameters, the
  selected Leiden or Louvain engine, graph names, cluster counts, and software
  versions in `@misc$isCORAL_metaclone`.
