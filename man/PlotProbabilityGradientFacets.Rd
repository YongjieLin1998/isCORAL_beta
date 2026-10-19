\name{PlotProbabilityGradientFacets}
\alias{PlotProbabilityGradientFacets}
\title{Plot pseudoclone probability gradients}
\usage{
PlotProbabilityGradientFacets(
  seurat_obj,
  demo_cells,
  highlight_colors,
  bg_gray = "#E0E4E7"
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing the pseudoclone cache.}
  \item{demo_cells}{Character vector of seed-cell names.}
  \item{highlight_colors}{Named colors for the displayed target lineages.}
  \item{bg_gray}{Color used for zero-probability cells.}
}
\value{A faceted ggplot object.}
\description{
Displays normalized pseudoclone probabilities for selected seed cells as
lineage-colored gradients across the cached cell embedding.
}

