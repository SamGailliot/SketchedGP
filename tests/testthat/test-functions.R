test_that("get_sketch_mat creates orthogonal matrix", {
  m = 3; p = 4
  M <- get_sketch_mat(3,4,orthog = TRUE)
  expect_equal(diag(M%*%t(M)), rep(1, m))
})

test_that("Correct kernel dimensions", {
  n1 <- 5; n2 <- 1; p = 10
  M1 <- get_sketch_mat(n1, p)
  M2 <- get_sketch_mat(n2, p)
  K1 <- exp_kernel(M1, lam = 1)
  K2 <- exp_kernel(M1, M2, lam = 1)
  K3 <- exp_kernel(M2, lam = 1)
  expect_equal(c(nrow(K1), ncol(K1)), c(n1, n1))
  expect_equal(c(nrow(K2), ncol(K2)), c(n1, n2))
  expect_equal(c(nrow(K3), ncol(K3)), c(n2, n2))
})
