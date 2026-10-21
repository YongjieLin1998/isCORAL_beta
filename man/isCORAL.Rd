\name{isCORAL}
\alias{isCORAL}
\title{Run the isCORAL beta 1.0.0 structural-correlation workflow}
\usage{
isCORAL(
  seurat_obj,
  assay = "RNA",
  input_mode = c("raw_counts", "pre_normalized_data"),
  pre_normalized_layer = "data",
  check_raw_counts = TRUE,
  results_dir = "iscoral_results_raw",
  species = "human",
  random_seed = 20251009,
  downsampleUMI = NULL,
  bin_size = 200,
  n_repeats = 100,
  parallel_cores = 25,
  max_cells = NULL,
  gemli_mean_q1 = 0.98,
  gemli_mean_q2 = 0.9,
  gemli_fano_q2 = 0.4,
  gemli_mean_q3 = 0.6,
  gemli_fano_q3 = 0.9,
  genes_to_remove = NULL,
  gene_min_mean_for_ranking = 0.05,
  nmf_feature_selection_ratio = 0.2,
  nmf_nrun = 1,
  nmf_feature_number = 5,
  nmf_k_value = 50,
  downsample_gemli_ratio = 1,
  specc_kernel = "laplacedot",
  nmf_eval_downsample_cells = NULL,
  custom_log_file = NULL,
  verbose = FALSE
)
}
\arguments{
  \item{seurat_obj}{A Seurat object.}
  \item{assay}{Assay containing the expression matrix.}
  \item{input_mode}{Input handling mode. The default, \code{"raw_counts"}, reads the counts layer and runs \code{Seurat::NormalizeData()}. Use \code{"pre_normalized_data"} only when the selected layer is already normalized.}
  \item{pre_normalized_layer}{Layer or slot to read when \code{input_mode = "pre_normalized_data"}.}
  \item{check_raw_counts}{Whether to warn when the counts layer appears processed rather than raw.}
  \item{results_dir}{Directory used for run logs and related output.}
  \item{species}{Species label used to identify mitochondrial and cell-cycle genes; currently \code{"human"} or \code{"mouse"}.}
  \item{random_seed}{Random seed used by the workflow and parallel backend.}
  \item{downsampleUMI}{Optional maximum UMI depth for raw-count input. Ignored for pre-normalized input.}
  \item{bin_size}{Target cell-bin size or vector of bin sizes.}
  \item{n_repeats}{Number of repeated structural scans per bin size.}
  \item{parallel_cores}{Number of CPU workers requested by the callr backend.}
  \item{max_cells}{Optional maximum number of cells retained before analysis.}
  \item{gemli_mean_q1, gemli_mean_q2, gemli_mean_q3}{Expression-mean quantiles used for GEMLI candidate-gene selection.}
  \item{gemli_fano_q2, gemli_fano_q3}{Fano-factor quantiles used for GEMLI candidate-gene selection.}
  \item{genes_to_remove}{Optional additional genes excluded from feature selection.}
  \item{gene_min_mean_for_ranking}{Minimum mean expression used in gene ranking.}
  \item{nmf_feature_selection_ratio}{Fraction of NMF features retained after ranking.}
  \item{nmf_nrun}{Number of NMF runs.}
  \item{nmf_feature_number}{Number of genes selected per NMF component.}
  \item{nmf_k_value}{Requested NMF rank.}
  \item{downsample_gemli_ratio}{Fraction of GEMLI genes sampled in each repeat.}
  \item{specc_kernel}{Kernel name passed to spectral clustering.}
  \item{nmf_eval_downsample_cells}{Optional number of cells used during NMF feature evaluation.}
  \item{custom_log_file}{Optional path for the progress log.}
  \item{verbose}{Whether to print routine progress and Seurat object-update messages. Defaults to \code{FALSE}. File logging, warnings, and errors remain enabled.}
}
\value{
A Seurat object with isCORAL scan results stored in
\code{seurat_obj@misc$CORAL_scan_results}. Run parameters record the package
version, input mode, normalization action, and expression layer used.
}
\description{
Runs the CPU isCORAL structural-correlation workflow. Beta 1.0.0 treats raw
counts as the default input, while allowing an explicit pre-normalized mode.
Assay access supports Seurat v4 slots and Seurat v5 layers.
}
\details{
The core scan uses a fork-based \code{doMC} backend and is intended for
Unix-like or HPC environments. A Windows-specific PSOCK backend is not
implemented in beta 1.0.0.
}
