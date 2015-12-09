.interp <- function(m, fits, slopes, at) {
  if(any(is.nan(fits))) {
    ind <- !is.nan(fits)
    c_interp(m[ind], fits[ind], slopes[ind], at)
  } else {
    c_interp(m, fits, slopes, at)
  }
}
