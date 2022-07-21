DISTMAT <- function (X){
  n <- ncol(X)
  m <- nrow(X)
  X <- as.matrix(X)
  A <- X %*% t(X)
  l <- diag(A)
  Dis <- sqrt(matrix(l, nrow = m, ncol =m)+ t(matrix(l, nrow = m, ncol =m)) - 2*A)
  return(Dis)
}
