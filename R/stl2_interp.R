.stl2_interp <- function(m, fits, slopes, at) {
   ind <- !is.nan(fits)
   out <- .C("stl2_interp", as.integer(m[ind]), as.double(fits[ind]), as.double(slopes[ind]), as.integer(length(fits[ind])), as.integer(at), as.integer(length(at)), ans=double(length(at)), PACKAGE="stl2")
   return(out$ans)
}
