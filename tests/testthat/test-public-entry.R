test_that("beta exposes only the unversioned core entry", {
  exports <- getNamespaceExports("isCORAL")

  expect_true("isCORAL" %in% exports)
  expect_false(any(grepl("^isCORAL_v[0-9_]+$", exports)))
  expect_true(exists("isCORAL", envir = asNamespace("isCORAL"), inherits = FALSE))
  expect_identical(as.character(utils::packageVersion("isCORAL")), "1.0.0")
})
