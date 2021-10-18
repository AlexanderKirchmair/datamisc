



test_that("ggvolcano", {

  testdf <- data.frame(row.names = as.vector(sapply(LETTERS, paste0, LETTERS)))
  testdf$lfc <- runif(nrow(testdf), min = -5, max = 5)
  testdf$fdr <- runif(nrow(testdf), min = 0, max = 1)

  gg <- ggvolcano(testdf)

  expect_true(all(class(gg) %in% c("gg", "ggplot")))
  expect_equal(nrow(gg$data), nrow(testdf))

})












