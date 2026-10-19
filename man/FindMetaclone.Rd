\name{FindMetaclone}
\alias{FindMetaclone}
\alias{FindMetaclones}
\title{Identify Metaclone communities in the CORAL clone space}
\usage{
FindMetaclone(
  seurat_obj,
  resolution = 0.7,
  reduction = "CORAL_clone_embed",
  dims = 1:30,
  graph_prefix = "pcoral",
  cluster_col_name = "metaclone",
  algorithm = 4,
  random_seed = 42,
  verbose = TRUE
)

FindMetaclones(
  seurat_obj,
  resolution = 0.7,
  reduction = "CORAL_clone_embed",
  dims = 1:30,
  graph_prefix = "pcoral",
  cluster_col_name = "metaclone",
  algorithm = 4,
  random_seed = 42,
  verbose = TRUE
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing the requested clone-space reduction.}
  \item{resolution}{Community-detection resolution passed to \code{Seurat::FindClusters()}.}
  \item{reduction}{Name of the CORAL clone-space reduction.}
  \item{dims}{Dimensions used to construct the nearest-neighbor graph. The
  historical default \code{1:30} is limited to the available dimensions for
  compatibility with lower-dimensional CORAL reductions; other out-of-range
  requests remain errors.}
  \item{graph_prefix}{Prefix for the generated nearest-neighbor and shared-nearest-neighbor graphs.}
  \item{cluster_col_name}{Metadata column receiving the Metaclone labels.}
  \item{algorithm}{Clustering algorithm passed to \code{Seurat::FindClusters()}; \code{4} requests Leiden and \code{1} requests Louvain.}
  \item{random_seed}{Random seed for community detection.}
  \item{verbose}{Whether to report graph construction and clustering progress.}
}
\value{The modified Seurat object with Metaclone labels, neighbor graphs, and
a run record in \code{@misc$isCORAL_metaclone}.}
\description{
Builds a shared-nearest-neighbor graph in the CORAL clone embedding and
identifies Metaclone communities with Leiden or Louvain clustering. Leiden
uses the optional \code{leidenbase} package through Seurat >= 5.2.0. If that
backend is unavailable, the function warns and falls back to Louvain.
}
\details{
\code{FindMetaclone()} is the canonical beta 1.0.0 entry point. The
deprecated \code{FindMetaclones()} name is retained only as a compatibility
wrapper for scripts written for versions 6.1.0 and 6.1.1. Both entry points
use the same implementation and parameter defaults. Optional software is
never installed while either function is running.
}
