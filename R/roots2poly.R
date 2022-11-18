roots2poly <- function(roots) {
  inv_roots <- 1 / roots

  n_coef <- length(roots)
  coefs <- numeric(n_coef)
  for (i in 1:n_coef) {
    coefs[i] <- utils::combn(x = inv_roots, m = i) |>
      apply(2, prod) |>
      sum() * (-1)^i
  }
  coefs
}


poly1 <- function(x) {
  val <- 1
  for (i in 1:length(roots)) {
    val <- val * (1 - inv_roots[i] * x)
  }
  val
}

poly2 <- function(x) {
  1 + my_poly[1] * x + my_poly[2] * x^2 + my_poly[3] * x^3
}


# 1 - sum(inv_roots) * x +
#   (inv_roots[1] * inv_roots[2] + inv_roots[1] * inv_roots[3] + inv_roots[2] * inv_roots[3]) * x^2 -
#   (prod(inv_roots)) * x^3
