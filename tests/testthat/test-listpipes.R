

test_that("List pipe operators", {

  A <- list("a" = 1:3, "b" = 1:5)

    expect_type(A %L>% sum(), "list")
  expect_equal(A %S>% sum(), c("a" = 6, "b" = 15))
  expect_equal(A %S>% function(x){ sum(x, 1:3) }, c("a" = 12, "b" = 21))
  expect_equal(A %S>% function(x){ .name }, c("a" = "a", "b" = "b"))


})


test_that("List assignment pipes", {

  A <- list("a" = 1:3, "b" = 1:5)
  Asqrt <- A %L>% sqrt()
  A %<L>% sqrt()

  expect_equal(Asqrt, A)

})


test_that("Multiple listpipes", {

  expect_equal((1:3 %L>% function(x){ c(x,x) }) %S>% sum(), c("1" = 2,"2" = 4,"3" = 6))


})


test_that("Parallel listpipes", {

  expect_equal(1:3 %L>% sqrt(), 1:3 %P>% sqrt())
  expect_equal(1:3 %L>% function(x){ .name }, 1:3 %P>% function(x){ .name })

})


test_that("Qsub pipe", {

  skip_if(.Platform$OS.type != "unix")
  skip_if(length(system(paste0("source ~/.bashrc; which qsub"), intern = TRUE)) == 0)
  expect_equal(1:3 %L>% sqrt(), 1:3 %Q>% sqrt())
  expect_equal(1:3 %L>% function(x){ .name }, 1:3 %Q>% function(x){ .name })

})



