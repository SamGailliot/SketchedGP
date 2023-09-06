test_that("get_sketch_mat creates orthogonal matrix", {
  m = 3; p = 4
  M <- get_sketch_mat(3,4,orthog = TRUE)
  expect_equal(diag(M%*%t(M)), rep(1, m))
})
