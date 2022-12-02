test_that("ARMApolyroots works", {

  set.seed(44322)
  inv_roots <- .sample_inv_roots(2)
  inv_roots <- c(inv_roots, Conj(inv_roots))
  ar_coef <- list(ar = Re(.roots2poly(inv_roots)))
  res <- 1 / ARMApolyroots(ar_coef)

  expect_equal(Re(res[order(Re(res))]), Re(inv_roots[order(Re(inv_roots))]))
  expect_equal(Im(res[order(Im(res))]), Im(inv_roots[order(Im(inv_roots))]))
})
