.ma <- function(x, n.p) {
  if(is.null(x) || !is.numeric(n.p))
    stop("something is wrong with the specified parameters")

  n <- length(x)
  out <- .C("ma", as.double(x), as.integer(n), as.integer(n.p),
    ans = double(n - 2 * n.p), PACKAGE = "stlplus")

  return(out$ans)
}
