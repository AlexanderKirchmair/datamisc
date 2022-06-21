


test_that("saveplot", {

  gg <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl)) + ggplot2::geom_point()
  file <- tempfile()
  min_fs <- 100

  saveplot(gg, file = file, dev = "png")
  on.exit(unlink(paste0(file, ".png")))
  expect_true(file.exists(paste0(file, ".png")))
  expect_gte(file.size(paste0(file, ".png")), min_fs)

  saveplot(gg, file = file, dev = "png", ggsave = FALSE)
  on.exit(unlink(paste0(file, ".png")))
  expect_true(file.exists(paste0(file, ".png")))
  expect_gte(file.size(paste0(file, ".png")), min_fs)

  saveplot(gg, file = file, dev = "tiff", ggsave = FALSE)
  on.exit(unlink(paste0(file, ".tiff")))
  expect_true(file.exists(paste0(file, ".tiff")))
  expect_gte(file.size(paste0(file, ".tiff")), min_fs)

  saveplot(gg, file = file, dev = "svg", ggsave = FALSE)
  on.exit(unlink(paste0(file, ".svg")))
  expect_true(file.exists(paste0(file, ".svg")))
  expect_gte(file.size(paste0(file, ".svg")), min_fs)

  saveplot(list(gg, gg), file = file, dev = "pdf", ggsave = FALSE)
  on.exit(unlink(paste0(file, ".pdf")))
  expect_true(file.exists(paste0(file, ".pdf")))
  expect_gte(file.size(paste0(file, ".pdf")), min_fs)

  saveplot(gg, file = file, dev = "cairo_ps")
  on.exit(unlink(paste0(file, ".cairo_ps")))
  expect_true(file.exists(paste0(file, ".cairo_ps")))
  expect_gte(file.size(paste0(file, ".cairo_ps")), min_fs)

})













