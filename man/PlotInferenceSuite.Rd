\name{PlotInferenceSuite}
\alias{PlotInferenceSuite}
\title{Plot the complete Metaclone inference matrix}
\usage{
PlotInferenceSuite(
  seurat_obj,
  reduction,
  space_title,
  cluster_col = "metaclone",
  custom_colors = NULL
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing clone-space coordinates, Metaclone labels, and cached pseudoclone probabilities.}
  \item{reduction}{Reduction used for the inference facets.}
  \item{space_title}{Text identifying the plotted space.}
  \item{cluster_col}{Metadata column containing Metaclone labels.}
  \item{custom_colors}{Optional named vector of cluster colors.}
}
\value{A patchwork object containing the complete inferred-probability facet matrix.}
\description{
Visualizes the transcriptome-inferred pseudoclone probability distribution of
every Metaclone, with seed cells shown in darker cluster colors.
}

