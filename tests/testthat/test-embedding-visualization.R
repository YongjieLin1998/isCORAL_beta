make_coral_embedding_plot_fixture <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("Matrix")
  cells <- paste0("cell_", 1:6)
  counts <- matrix(
    2, nrow = 3, ncol = 6,
    dimnames = list(paste0("gene_", 1:3), cells)
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  Seurat::Idents(object) <- factor(c("A", "B", "A", "B", "A", "B"))
  object$cell_type <- factor(c("B", "A", "B", NA, "C", "A"),
                             levels = c("C", "A", "B"))
  object$metaclone <- c("2", "2", "3", "3", "4", "4")
  coordinates <- cbind(
    CORALcell_1 = 1:6, CORALcell_2 = c(1, 3, 2, 5, 4, 6),
    CORALcell_3 = c(6, 5, 4, 3, 2, 1)
  )
  rownames(coordinates) <- cells
  object[["CORAL_cell_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = coordinates, key = "CORALcell_", assay = "RNA"
  )
  clone_coordinates <- coordinates / 2
  colnames(clone_coordinates) <- paste0("CORALclone_", 1:3)
  object[["CORAL_clone_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = clone_coordinates, key = "CORALclone_", assay = "RNA"
  )
  probabilities <- diag(6)
  probabilities[1, ] <- c(0.4, 0.2, 0.3, 0.1, 0, 0)
  probabilities[2, ] <- c(0.1, 0.5, 0.2, 0.1, 0.1, 0)
  dimnames(probabilities) <- list(cells, cells)
  object@misc$isCORAL_embedding <- list(
    P_norm = Matrix::Matrix(probabilities, sparse = TRUE),
    valid_cells = cells,
    source_binsize = "binsize_200"
  )
  object
}

test_that("PlotCORALEmbedding preserves coordinates, identities, and the input", {
  object <- make_coral_embedding_plot_fixture()
  original <- object
  coordinates <- Seurat::Embeddings(object, "CORAL_cell_embed")
  p <- PlotCORALEmbedding(object)

  expect_s3_class(p, "ggplot")
  expect_identical(p$data$Cell_ID, rownames(coordinates))
  expect_equal(p$data$Dim1, unname(coordinates[, 1]))
  expect_equal(p$data$Dim2, unname(coordinates[, 2]))
  expect_equal(as.character(p$data$Group), as.character(Seurat::Idents(object)))
  expect_equal(p$layers[[1]]$aes_params$shape, 16)
  expect_equal(p$theme$text$family, "Arial")
  expect_equal(p$theme$text$size, 18)
  expect_identical(object, original)
})

test_that("PlotCORALEmbedding defaults to hollow clone points with an override", {
  object <- make_coral_embedding_plot_fixture()
  p <- PlotCORALEmbedding(object, reduction = "CORAL_clone_embed")
  filled <- PlotCORALEmbedding(
    object, reduction = "CORAL_clone_embed", hollow = FALSE
  )
  expect_equal(p$layers[[1]]$aes_params$shape, 1)
  expect_equal(filled$layers[[1]]$aes_params$shape, 16)
  expect_equal(p$data, filled$data)
})

test_that("PlotCORALEmbedding respects metadata factors, missing labels, and dimensions", {
  object <- make_coral_embedding_plot_fixture()
  coordinates <- Seurat::Embeddings(object, "CORAL_clone_embed")
  p <- PlotCORALEmbedding(
    object, reduction = "CORAL_clone_embed", group_by = "cell_type",
    dims = c(3, 1), colors = c(C = "#3FA989", A = "#D88DBA", B = "#F0AE46")
  )
  expect_identical(levels(p$data$Group), c("C", "A", "B"))
  expect_true(is.na(p$data$Group[4]))
  expect_equal(nrow(p$data), ncol(object))
  expect_equal(p$data$Dim1, unname(coordinates[, 3]))
  expect_equal(p$data$Dim2, unname(coordinates[, 1]))
  expect_equal(p$scales$get_scales("colour")$na.value, "#A9A9A9")
})

test_that("PlotCORALEmbedding aligns subsetted reduction cells by name", {
  object <- make_coral_embedding_plot_fixture()
  cells <- colnames(object)[c(2, 4, 6)]
  coordinates <- Seurat::Embeddings(object, "CORAL_cell_embed")[cells, , drop = FALSE]
  colnames(coordinates) <- paste0("Subset_", 1:3)
  object[["subset_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = coordinates, key = "Subset_", assay = "RNA"
  )
  p <- PlotCORALEmbedding(object, reduction = "subset_embed", group_by = "cell_type")
  expect_identical(p$data$Cell_ID, cells)
  expect_equal(as.character(p$data$Group), c("A", NA, "A"))
})

test_that("PlotCORALEmbedding reports invalid plotting inputs", {
  object <- make_coral_embedding_plot_fixture()
  expect_error(PlotCORALEmbedding(NULL), "Seurat object")
  expect_error(PlotCORALEmbedding(object, reduction = "missing"), "reduction was not found")
  expect_error(PlotCORALEmbedding(object, group_by = "missing"), "metadata column")
  expect_error(PlotCORALEmbedding(object, dims = c(1, 1)), "two different")
  expect_error(PlotCORALEmbedding(object, dims = c(1, 4)), "two different")
  expect_error(PlotCORALEmbedding(object, dims = c(1, Inf)), "two different")
  expect_error(PlotCORALEmbedding(object, alpha = 2), "between 0 and 1")
  expect_error(PlotCORALEmbedding(object, hollow = NA), "TRUE, or FALSE")
  expect_error(PlotCORALEmbedding(object, point_size = 0), "positive number")
  expect_error(PlotCORALEmbedding(object, colors = c(A = "red")), "every observed group")
})

test_that("Pseudoclone clone-space circles retain the same enrichment calculation", {
  object <- make_coral_embedding_plot_fixture()
  p_cell <- PlotPseudoclone(
    object, seed_cell = "cell_1", reduction = "CORAL_cell_embed",
    prob_mat = object@misc$isCORAL_embedding$P_norm,
    min_enrichment = 0.5, show_ground_truth = FALSE, return_data = TRUE
  )
  p_clone <- PlotPseudoclone(
    object, seed_cell = "cell_1", reduction = "CORAL_clone_embed",
    prob_mat = object@misc$isCORAL_embedding$P_norm,
    min_enrichment = 0.5, show_ground_truth = FALSE, return_data = TRUE
  )
  expect_equal(p_clone$prediction_data$Enrichment, p_cell$prediction_data$Enrichment)
  expect_equal(p_clone$prediction_data$Assigned_Lineage, p_cell$prediction_data$Assigned_Lineage)
  expect_equal(p_clone$seed_cells, p_cell$seed_cells)
  expect_equal(p_clone$plot$layers[[1]]$aes_params$shape, 1)
  expect_equal(p_clone$plot$layers[[2]]$aes_params$shape, 1)
  expect_equal(p_clone$plot$layers[[3]]$aes_params$shape, 21)
  expect_equal(p_cell$plot$layers[[1]]$aes_params$shape, 19)
  expect_equal(p_cell$plot$layers[[2]]$aes_params$shape, 19)
})

test_that("Metaclone clone-space circles preserve probabilities and filled seeds", {
  object <- make_coral_embedding_plot_fixture()
  p_cell <- PlotMetaclone(object, 2, return_data = TRUE)
  p_clone <- PlotMetaclone(
    object, 2, reduction = "CORAL_clone_embed", return_data = TRUE
  )
  expect_equal(p_clone$plot_data$RelativeProbability, p_cell$plot_data$RelativeProbability)
  expect_equal(p_clone$seed_cells, p_cell$seed_cells)
  expect_equal(p_clone$sister_cells, p_cell$sister_cells)
  expect_equal(p_clone$plot$layers[[1]]$aes_params$shape, 1)
  expect_equal(p_clone$plot$layers[[2]]$aes_params$shape, 1)
  expect_equal(p_clone$plot$layers[[3]]$aes_params$shape, 16)
  expect_equal(p_cell$plot$layers[[1]]$aes_params$shape, 16)
  expect_equal(p_cell$plot$layers[[2]]$aes_params$shape, 16)
})
