\name{PlotValidationSuite}
\alias{PlotValidationSuite}
\title{Plot inferred and experimental clone validation facets}
\usage{
PlotValidationSuite(
  seurat_obj,
  reduction,
  space_title,
  barcode_col,
  cluster_col = "metaclone",
  na_barcodes = NULL,
  custom_colors = NULL
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing clone-space coordinates, Metaclone labels, and cached pseudoclone probabilities.}
  \item{reduction}{Reduction used for the validation panels.}
  \item{space_title}{Text identifying the plotted space.}
  \item{barcode_col}{Metadata column containing experimental clone labels.}
  \item{cluster_col}{Metadata column containing Metaclone labels.}
  \item{na_barcodes}{Optional additional barcode labels treated as invalid or background.}
  \item{custom_colors}{Optional named vector of cluster colors.}
}
\value{A patchwork object comparing inferred pseudoclone structure with experimental clone structure.}
\description{
Builds matched faceted panels for transcriptome-inferred pseudoclone
probabilities and experimental lineage labels. Dense point layers are
rasterized for efficient rendering.
}

