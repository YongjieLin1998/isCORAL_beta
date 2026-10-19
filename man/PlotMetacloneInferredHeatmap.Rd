\name{PlotMetacloneInferredHeatmap}
\alias{PlotMetacloneInferredHeatmap}
\title{Plot inferred Metaclone lineage potential}
\usage{
PlotMetacloneInferredHeatmap(
  seurat_obj,
  cluster_col = "metaclone",
  cell_type_col = "cluster_integrated",
  target_lineages = NULL
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing Metaclone labels, cell types, and cached pseudoclone probabilities.}
  \item{cluster_col}{Metadata column containing Metaclone labels.}
  \item{cell_type_col}{Metadata column containing cell-type labels.}
  \item{target_lineages}{Optional ordered character vector of lineages to display.}
}
\value{A list with the inferred percentage \code{table} and combined \code{plot}.}
\description{
Aggregates cached pseudoclone probabilities by Metaclone and target lineage,
then displays the inferred fate matrix with a clustering tree.
}

