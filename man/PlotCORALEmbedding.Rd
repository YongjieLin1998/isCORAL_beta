\name{PlotCORALEmbedding}
\alias{PlotCORALEmbedding}
\title{Plot a pCORAL cell or pseudoclone embedding}
\usage{
PlotCORALEmbedding(
  seurat_obj,
  reduction = "CORAL_cell_embed",
  group_by = NULL,
  dims = c(1, 2),
  colors = NULL,
  point_size = 1.2,
  alpha = 1,
  hollow = NULL,
  stroke = 0.6,
  base_size = 18,
  base_family = "Arial",
  title = NULL
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing the requested reduction.}
  \item{reduction}{Reduction to plot. Defaults to \code{"CORAL_cell_embed"}; use \code{"CORAL_clone_embed"} for the pseudoclone embedding.}
  \item{group_by}{Optional metadata column for categorical coloring. \code{NULL} uses the current Seurat identities.}
  \item{dims}{Two embedding dimensions to display.}
  \item{colors}{Optional color vector, preferably named by group.}
  \item{point_size}{Point size.}
  \item{alpha}{Point opacity, between 0 and 1.}
  \item{hollow}{Use open circles. \code{NULL} selects open circles for \code{"CORAL_clone_embed"} and filled points for other reductions.}
  \item{stroke}{Outline width for open circles.}
  \item{base_size, base_family}{Base font size and font family.}
  \item{title}{Optional plot title.}
}
\description{
Displays existing embedding coordinates with Arial text and no axes.
Cell embeddings use filled points; clone embeddings use open circles.
}
\details{
Coordinates and grouping labels are aligned by cell name. All cells in the
selected reduction are plotted without subsampling. Numeric metadata are
treated as categorical labels, and missing labels remain missing and are
shown in gray. Named colors must cover all observed non-missing groups.

This function reads only the reduction and cell labels. It does not
recompute an embedding, access expression layers, modify the Seurat object,
print a figure, or save a file.
}
\value{
A ggplot object. Its \code{data} contains \code{Cell_ID}, \code{Dim1},
\code{Dim2}, and \code{Group}.
}
\examples{
\dontrun{
p_cell <- PlotCORALEmbedding(seurat_obj)
print(p_cell)

p_clone <- PlotCORALEmbedding(seurat_obj, reduction = "CORAL_clone_embed")
print(p_clone)

p_metaclone <- PlotCORALEmbedding(
  seurat_obj, reduction = "CORAL_clone_embed", group_by = "metaclone"
)
print(p_metaclone)
}
}
