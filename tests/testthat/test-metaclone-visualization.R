make_metaclone_plot_fixture <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("Matrix")
  cells <- paste0("cell_", 1:6)
  counts <- matrix(2, nrow = 3, ncol = 6,
                   dimnames = list(paste0("gene_", 1:3), cells))
  object <- Seurat::CreateSeuratObject(counts = counts)
  object$metaclone <- c("2", "2", "3", "3", "4", "4")
  object$cell_type <- rep("TypeA", 6)
  coordinates <- cbind(CORALcell_1 = 1:6, CORALcell_2 = c(1, 3, 2, 5, 4, 6))
  rownames(coordinates) <- cells
  object[["CORAL_cell_embed"]] <- Seurat::CreateDimReducObject(
    embeddings = coordinates, key = "CORALcell_", assay = "RNA"
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

test_that("PlotMetaclone preserves the embedding script's probability calculation", {
  object <- make_metaclone_plot_fixture()
  original_cache <- object@misc$isCORAL_embedding
  result <- PlotMetaclone(object, metaclone_id = "2", return_data = TRUE)

  expect_s3_class(result$plot, "ggplot")
  expect_equal(result$plot_data$MeanProbability, c(0.25, 0.35, 0.25, 0.1, 0.05, 0))
  expect_equal(result$plot_data$SisterProbability, c(0, 0, 0.25, 0.1, 0.05, 0))
  expect_equal(result$plot_data$RelativeProbability, c(0, 0, 1, 0.4, 0.2, 0))
  expect_identical(result$seed_cells, c("cell_1", "cell_2"))
  expect_identical(result$sister_cells, c("cell_3", "cell_4", "cell_5"))
  expect_identical(result$reduction, "CORAL_cell_embed")
  expect_identical(result$source_binsize, "binsize_200")
  expect_identical(object@misc$isCORAL_embedding, original_cache)
})

test_that("PlotMetaclone retains the gray underlay, ordered cloud, and dark seeds", {
  object <- make_metaclone_plot_fixture()
  result <- PlotMetaclone(object, metaclone_id = 2, return_data = TRUE)
  layers <- result$plot$layers

  expect_length(layers, 3)
  expect_identical(layers[[1]]$data$Cell_ID, colnames(object))
  expect_identical(layers[[2]]$data$Cell_ID, c("cell_5", "cell_4", "cell_3"))
  expect_identical(layers[[3]]$data$Cell_ID, c("cell_1", "cell_2"))
  expect_equal(layers[[1]]$aes_params$colour, "#E0E4E7")
  expect_equal(layers[[2]]$aes_params$colour, "#3182BD")
  expect_equal(layers[[3]]$aes_params$colour, result$plot_parameters$seed_color)
  expect_equal(layers[[3]]$aes_params$size, 2)
  expect_equal(layers[[3]]$aes_params$shape, 16)
})

test_that("PlotMetaclone aligns independently reordered matrix axes by cell name", {
  object <- make_metaclone_plot_fixture()
  expected <- PlotMetaclone(object, 2, return_data = TRUE)$plot_data
  p <- object@misc$isCORAL_embedding$P_norm
  object@misc$isCORAL_embedding$P_norm <- p[c(4, 2, 6, 1, 3, 5), c(6, 5, 3, 1, 2, 4)]
  actual <- PlotMetaclone(object, 2, return_data = TRUE)$plot_data
  expect_equal(actual, expected)
})

test_that("PlotMetaclone uses valid_cells for unnamed cache axes", {
  object <- make_metaclone_plot_fixture()
  expected <- PlotMetaclone(object, 2, return_data = TRUE)$plot_data
  dimnames(object@misc$isCORAL_embedding$P_norm) <- list(NULL, NULL)
  actual <- PlotMetaclone(object, 2, return_data = TRUE)$plot_data
  expect_equal(actual, expected)
})

test_that("PlotMetaclone accepts character labels and custom metadata columns", {
  object <- make_metaclone_plot_fixture()
  object$custom_metaclone <- c("Stem", "Stem", "Other", "Other", "Other", "Other")
  result <- PlotMetaclone(object, "Stem", cluster_col = "custom_metaclone", return_data = TRUE)
  expect_identical(result$seed_cells, c("cell_1", "cell_2"))
  expect_identical(result$metaclone_id, "Stem")
})

test_that("PlotMetaclone handles a zero sister field without replacing it by NA", {
  object <- make_metaclone_plot_fixture()
  p <- diag(6)
  dimnames(p) <- list(colnames(object), colnames(object))
  object@misc$isCORAL_embedding$P_norm <- p
  result <- PlotMetaclone(object, 2, return_data = TRUE)
  expect_equal(result$plot_data$RelativeProbability, rep(0, 6))
  expect_length(result$sister_cells, 0)
  expect_equal(nrow(result$plot$layers[[2]]$data), 0)
})

test_that("PlotMetaclone preserves unscored cells as background with NA values", {
  object <- make_metaclone_plot_fixture()
  object@misc$isCORAL_embedding$P_norm <- object@misc$isCORAL_embedding$P_norm[1:5, 1:5]
  object@misc$isCORAL_embedding$valid_cells <- colnames(object)[1:5]
  expect_warning(
    result <- PlotMetaclone(object, 2, return_data = TRUE),
    "no cached probability"
  )
  expect_true(is.na(result$plot_data$RelativeProbability[6]))
  expect_false(result$plot_data$IsSister[6])
  expect_equal(nrow(result$plot$layers[[1]]$data), 6)
})

test_that("PlotMetaclone adds the optional embedding-style contour layer", {
  object <- make_metaclone_plot_fixture()
  result <- PlotMetaclone(
    object, 2, cell_type_col = "cell_type",
    contour_colors = c(TypeA = "#D88DBA"), return_data = TRUE
  )
  expect_length(result$plot$layers, 4)
  expect_equal(result$plot$layers[[4]]$aes_params$linetype, "dashed")
  expect_equal(result$plot_parameters$contour_level, 0.15)
})

test_that("PlotMetaclone reports missing seeds and invalid inputs", {
  object <- make_metaclone_plot_fixture()
  expect_error(PlotMetaclone(object, "missing"), "was not found")
  expect_error(PlotMetaclone(object, c(2, 3)), "one non-missing")
  expect_error(PlotMetaclone(object, 2, dims = c(1, 1)), "two different")
  expect_error(PlotMetaclone(object, 2, reduction = "CORAL_cell_embedding"), "reduction was not found")
  expect_error(PlotMetaclone(object, 2, cell_type_col = "cell_type"), "contour_colors")
  object@misc$isCORAL_embedding$P_norm <- object@misc$isCORAL_embedding$P_norm[-2, -2]
  expect_error(PlotMetaclone(object, 2), "seed cell.*missing")
})

test_that("PlotMetaclone rejects missing or non-finite probability caches", {
  object <- make_metaclone_plot_fixture()
  object@misc$isCORAL_embedding$P_norm[1, 3] <- NA_real_
  expect_error(PlotMetaclone(object, 2), "non-finite or negative")
  object@misc$isCORAL_embedding$P_norm <- NULL
  expect_error(PlotMetaclone(object, 2), "run CORALembedding")
})
