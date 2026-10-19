\name{PlotProbabilityCollapseDiptych}
\alias{PlotProbabilityCollapseDiptych}
\title{Compare inferred probabilities with experimental siblings}
\usage{
PlotProbabilityCollapseDiptych(
  seurat_obj,
  demo_cells,
  highlight_colors,
  clone_col,
  bg_gray = "#E0E4E7",
  point_size = 0.5
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing the pseudoclone cache and a \code{umap.raw} reduction.}
  \item{demo_cells}{Character vector of seed-cell names.}
  \item{highlight_colors}{Named colors for the displayed cell lineages.}
  \item{clone_col}{Metadata column containing experimental clone labels.}
  \item{bg_gray}{Color used for background cells.}
  \item{point_size}{Base point size.}
}
\value{A vertically arranged patchwork object with one inferred/experimental diptych per valid seed cell.}
\description{
For selected seed cells, contrasts isCORAL pseudoclone probabilities with
experimentally labeled sibling cells.
}

