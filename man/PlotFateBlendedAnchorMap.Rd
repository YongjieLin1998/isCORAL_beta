\name{PlotFateBlendedAnchorMap}
\alias{PlotFateBlendedAnchorMap}
\title{Plot fate-blended pseudoclone anchors}
\usage{
PlotFateBlendedAnchorMap(
  seurat_obj,
  demo_cells,
  highlight_colors,
  reduction = "CORAL_clone_embed"
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing the requested reduction and pseudoclone cache.}
  \item{demo_cells}{Character vector of seed-cell names.}
  \item{highlight_colors}{Named colors for HSC, GMP, Mk, and Other lineages.}
  \item{reduction}{Reduction used to display clone-space coordinates.}
}
\value{A ggplot object.}
\description{
Maps selected seed cells onto the clone space and blends lineage colors
according to their inferred pseudoclone probability distributions.
}

