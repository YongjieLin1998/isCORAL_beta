\name{PlotMetacloneLandscape}
\alias{PlotMetacloneLandscape}
\title{Plot the Metaclone landscape}
\usage{
PlotMetacloneLandscape(
  seurat_obj,
  reduction = "CORAL_clone_embed",
  cluster_col = "metaclone",
  custom_colors = NULL
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing the requested reduction and Metaclone labels.}
  \item{reduction}{Reduction used for the two-dimensional landscape.}
  \item{cluster_col}{Metadata column containing Metaclone labels.}
  \item{custom_colors}{Optional named vector of cluster colors.}
}
\value{A patchwork object containing the landscape and cluster facets.}
\description{
Draws the Metaclone embedding, density contours, and per-community facets with
a consistent color palette.
}

