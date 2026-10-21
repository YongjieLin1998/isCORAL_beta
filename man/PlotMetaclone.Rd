\name{PlotMetaclone}
\alias{PlotMetaclone}
\title{Plot the seeds and inferred sisters of one Metaclone}
\usage{
PlotMetaclone(
  seurat_obj,
  metaclone_id,
  reduction = "CORAL_cell_embed",
  cluster_col = "metaclone",
  dims = c(1, 2),
  color = "#3182BD",
  seed_darken = 0.5,
  background_color = "#E0E4E7",
  background_size = 0.8,
  point_size = 1.0,
  seed_size = 2.0,
  cell_type_col = NULL,
  contour_colors = NULL,
  contour_level = 0.15,
  base_size = 18,
  base_family = "Arial",
  title = NULL,
  return_data = FALSE
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing Metaclone metadata and the probability cache produced by \code{CORALembedding()}.}
  \item{metaclone_id}{One Metaclone label from \code{cluster_col}.}
  \item{reduction}{Reduction to plot. Defaults to \code{"CORAL_cell_embed"}. In \code{"CORAL_clone_embed"}, background and sister cells use open circles; seeds remain filled.}
  \item{cluster_col}{Metadata column containing Metaclone labels.}
  \item{dims}{Two embedding dimensions to display.}
  \item{color}{Color of the inferred sister-cell cloud.}
  \item{seed_darken}{Multiplier applied to the cloud color's RGB channels to obtain the seed color. Use 0.5 for the single-Metaclone embedding example or 0.45 for the two-Metaclone example.}
  \item{background_color}{Color of the complete background cell layer.}
  \item{background_size, point_size, seed_size}{Point sizes for background, inferred sisters, and seeds.}
  \item{cell_type_col}{Optional metadata column used for density contours. Supply together with \code{contour_colors}.}
  \item{contour_colors}{Named colors for the cell types to contour. Each selected type needs at least three plotted cells and variation on both axes.}
  \item{contour_level}{Normalized density contour level.}
  \item{base_size, base_family}{Base font size and font family.}
  \item{title}{Plot title. \code{NULL} uses the selected Metaclone label.}
  \item{return_data}{If \code{TRUE}, return the plot, cell-level probabilities, seed and sister cell IDs, and plotting parameters.}
}
\description{
Uses the embedding scripts' probability-cloud style: a complete gray cell
background, a colored inferred sister cloud with probability-dependent
opacity, and darker seed cells.
}
\details{
All cells assigned to the selected Metaclone are seeds. Their rows of
\code{seurat_obj@misc$isCORAL_embedding$P_norm} are averaged, the probabilities
at seed cells are set to zero, and the result is divided by its maximum
remaining value. All non-seed cells with positive values contribute to the
inferred sister cloud. This preserves the calculation in the original
\code{embedding/5.6.R} and \code{embedding/5.12.R} scripts.

Metadata, matrix axes, and coordinates are aligned by cell name. An unnamed
cache axis uses \code{cache$valid_cells}. Averaging and maximum scaling use
all cached target cells, including when the object has been subsetted.
Cells with coordinates but no cached probability remain gray and retain
\code{NA} probabilities in the returned table. Optional contours use dashed
normalized-density lines, with a default level of 0.15.

The function uses the existing embedding and does not modify the Seurat
object, print a figure, or save a file.
}
\value{
A ggplot object. With \code{return_data = TRUE}, a list containing
\code{plot}, \code{plot_data}, \code{seed_cells}, \code{sister_cells},
\code{metaclone_id}, \code{reduction}, \code{cluster_col},
\code{source_binsize}, and \code{plot_parameters}.
The table includes \code{MeanProbability} before seed exclusion,
\code{SisterProbability} after seed exclusion, and
\code{RelativeProbability} after maximum scaling.
}
\examples{
\dontrun{
p <- PlotMetaclone(seurat_obj, metaclone_id = "2")
print(p)

p <- PlotMetaclone(
  seurat_obj, metaclone_id = "2",
  cell_type_col = "cluster_integrated",
  contour_colors = c(HSC = "#D88DBA", Mk = "#F0AE46", GMP = "#3FA989")
)
print(p)
}
}
