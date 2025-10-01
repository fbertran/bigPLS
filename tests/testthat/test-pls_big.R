test_that("pls_big handles in-memory big.matrix", {
  set.seed(123)
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  Y <- matrix(rnorm(40), nrow = 20, ncol = 2)
  bm <- bigmemory::as.big.matrix(X)
  res <- pls_big(bm, Y, ncomp = 2)
  expect_named(res, c("scores", "Yscores", "weights", "loadings", "Yloadings", "coefficients"))
  expect_equal(dim(res$scores), c(20, 2))
  expect_equal(dim(res$weights), c(3, 2))
  expect_equal(length(res$coefficients), 2)
})

test_that("streaming PLS matches compiled implementation for file-backed matrices", {
  skip_on_cran()
  set.seed(42)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  Y <- matrix(rnorm(60), nrow = 30, ncol = 2)
  tmpdir <- tempdir()
  backingfile <- "pls_stream.bin"
  descriptorfile <- "pls_stream.desc"
  fbm <- bigmemory::filebacked.big.matrix(
    nrow = nrow(X), ncol = ncol(X), type = "double",
    backingpath = tmpdir, backingfile = backingfile,
    descriptorfile = descriptorfile, init = X
  )
  bigmemory::flush(fbm)
  desc <- dget(file.path(tmpdir, descriptorfile))
  on.exit({
    try(unlink(file.path(tmpdir, backingfile)))
    try(unlink(file.path(tmpdir, descriptorfile)))
  }, add = TRUE)
  attached <- bigmemory::attach.big.matrix(desc)
  res_cpp <- pls_big(attached, Y, ncomp = 2)
  res_stream <- matrixpls_stream_bigmatrix(attached, Y, ncomp = 2, num.rows.chunk = 10)
  expect_equal(res_cpp$coefficients, res_stream$coefficients, tolerance = 1e-6)
  expect_equal(res_cpp$weights, res_stream$weights, tolerance = 1e-5)
  expect_equal(res_cpp$loadings, res_stream$loadings, tolerance = 1e-5)
})
