\name{PlotMetacloneSeedHeatmap}
\alias{PlotMetacloneSeedHeatmap}
\title{Plot Metaclone seed-cell composition}
\usage{
PlotMetacloneSeedHeatmap(
  seurat_obj,
  cluster_col = "metaclone",
  cell_type_col = "cluster_integrated",
  target_lineages = NULL
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing Metaclone and cell-type metadata.}
  \item{cluster_col}{Metadata column containing Metaclone labels.}
  \item{cell_type_col}{Metadata column containing cell-type labels.}
  \item{target_lineages}{Optional ordered character vector of lineages to display.}
}
\value{A list with the clustered count \code{table} and combined \code{plot}.}
\description{
Summarizes the observed seed-cell composition of each Metaclone and displays
the result with a hierarchical clustering tree.
}

