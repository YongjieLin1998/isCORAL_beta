\name{PlotMetacloneGroundTruthHeatmap}
\alias{PlotMetacloneGroundTruthHeatmap}
\title{Plot experimental sibling-lineage ground truth}
\usage{
PlotMetacloneGroundTruthHeatmap(
  seurat_obj,
  cluster_col = "metaclone",
  cell_type_col = "cluster_integrated",
  barcode_col,
  target_lineages = NULL,
  na_barcodes = c("NA", "Unknown_NA", "6", NA, "")
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing Metaclone, cell-type, and experimental clone metadata.}
  \item{cluster_col}{Metadata column containing Metaclone labels.}
  \item{cell_type_col}{Metadata column containing cell-type labels.}
  \item{barcode_col}{Metadata column containing experimental clone labels.}
  \item{target_lineages}{Optional ordered character vector of lineages to display.}
  \item{na_barcodes}{Barcode labels treated as invalid or background.}
}
\value{A list with the sibling-fate percentage \code{table} and combined \code{plot}.}
\description{
Calculates the experimental sibling-lineage distribution associated with each
Metaclone and displays the resulting ground-truth matrix with a clustering
tree.
}

