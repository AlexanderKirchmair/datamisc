

test_that("head2", {
  m <- rmat(100, 100)
  expect_equal( head2(m, nrows = 3, ncols = 5), m[1:3,1:5] )
})




test_that("NA handling", {

  mat <- matrix(1:12, 3) > 5
  mat[c(1,4,7)] <- NA

  expect_equal(sum(is.na(naf(mat))), 0)
  expect_equal(sum(is.na(nat(mat))), 0)

  expect_equal(sum(naf(mat)), 6)
  expect_equal(sum(nat(mat)), 9)

  expect_equal(naSkip(mat), mat[!is.na(rowSums(mat)),])

})



test_that("rmat", {
  expect_equal( dim(rmat(2,2)), c(2,2) )
  expect_true( all(rmat(FUN = rbinom, size = 1, prob = 1) == 1))
})


test_that("untidy rownames", {
  m0 <- rmat()
  m1 <- rownames2col(m0)
  m2 <- col2rownames(m1)
  expect_equal(as.matrix(m0), as.matrix(m2))
})


test_that("cjoin", {
  expect_equal( ncol(cjoin(mtcars, mtcars, mtcars)), ncol(mtcars)*3 )
  expect_equal( nrow(cjoin(mtcars[1:10,], mtcars)), 10 )
  expect_equal( nrow(cjoin(mtcars[1:10,], mtcars[6:15,])), 5)
  expect_equal(colnames(cjoin(mtcars[,c("mpg", "cyl")], NEW = setNames(mtcars[,1], rownames(mtcars)))), c("mpg", "cyl", "NEW"))
})


test_that("rjoin", {
  expect_equal( rjoin(mtcars[1:3,], mtcars[4:10,], mtcars[11:32,]), mtcars )
  expect_equal( dim(rjoin(mtcars[1:3,1:5], mtcars[4:10,3:8])), c(10, 3))
})


test_that("matScale", {

  expect_equal(dimnames(matScale(mtcars)), dimnames(mtcars))

  df <- data.frame(t(apply(iris[,-5], 1, scale)))
  dimnames(df) <- dimnames(iris[,-5])
  expect_equal(matScale(iris, rows = TRUE)[,-5], df)

  df <- data.frame(apply(iris[,-5], 2, scale))
  dimnames(df) <- dimnames(iris[,-5])
  expect_equal(matScale(iris, cols = TRUE)[,-5], df)

})


test_that("stratify", {
  expect_equal(as.numeric(stratify(1:10)), ifelse(1:10 <= median(1:10), 1, 2))
  expect_equal(as.numeric(stratify(1:6, 3)), c(1,1,2,2,3,3))
  expect_equal(as.numeric(stratify(mtcars$cyl)), ifelse(mtcars$cyl <= median(mtcars$cyl), 1, 2))
})


test_that("dedupl", {
  x <- paste0("aa", LETTERS[c(1,1,1,2,3,4,5,5)])
  expect_equal(sum(duplicated(dedupl(x))), 0)
  expect_equal(sum(grepl("__", dedupl(x, sep = "__", index = LETTERS))), sum(duplicated(x)))
})


test_that("padjust", {
  pvec <- log10(runif(100, 1, 10))
  pvec[sample(seq(pvec), size = 10)] <- 0.00001
  expect_equal(padjust(pvec, method = "holm"), p.adjust(pvec, method = "holm"))
  expect_equal(padjust(pvec, method = "fdr"), p.adjust(pvec, method = "fdr"))
  expect_equal(as.vector(padjust(matrix(pvec, nrow = 10), method = "fdr")), p.adjust(pvec, method = "fdr"))
})



test_that("Write and read tables", {
  filename <- "testfile_testthat.xlsx"
  input <- list("mtcars" = mtcars, "iris" = iris)
  writeTables(input, file = filename)
  output <- readTables(file = filename)
  unlink(filename)

  expect_equal(names(output), names(input))
  expect_equal(dim(output$mtcars), dim(input$mtcars))
  expect_equal(colnames(output$mtcars), colnames(input$mtcars))
  expect_equal(rownames(output$mtcars), rownames(input$mtcars))

  input <- list("iris" = iris, "iris" = iris)
  writeTables(input, file = filename)
  output <- readTables(file = filename)
  unlink(filename)

  expect_equal(sum(duplicated(names(input))), 1)
  expect_equal(sum(duplicated(names(output))), 0)
})



