#' Accessor functions for elements of an stl and stlplus object
#'
#'Retrieves the raw, seasonal, trend, remainder, or time components from an stlplus object.  The methods \code{seasonal.stl}, ... also exist as a convenience for extracting components from R's \code{stl()}.
#'
#' @param x,object object of class \code{"stl"} or \code{"stlplus"}.
#' @param fcnum number of post-trend smoothing frequency component.
#' @param \ldots additional parameters
#'
#' @return
#' Returns a vector of either the \code{getraw} time series, the \code{seasonal}, \code{trend}, or \code{remainder} components, or the \code{time} values of the time series.  If \code{time}s are requested but were not supplied in the initial \code{stlplus} call, the \code{1:n} vector is returned, where \code{n} is the number of data points.  The \code{fitted} method returns the sum of the seasonal and trend.
#' @references R. B. Cleveland, W. S. Cleveland, J.E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
#' @note
#' The \code{fitted} and \code{predict} methods are equivalent.  For objects of class \code{"stlplus"}, these functions return the sum of all components but the remainder, including post-trend smoothing components.  Note also that the \code{trend} method for objects of class \code{"stlplus"} only returns the trend component from the STL iterations, even when post-trend smoothing is done.
#' @seealso \code{\link{stlplus}}
#' @examples
#' co2.stl <- stlplus(co2, t = as.vector(stats::time(co2)), n.p=12, l.window=13,
#' t.window=19, s.window=35, s.degree=1, sub.labels = substr(month.name, 1, 3))
#'
#' plot(seasonal(co2.stl))
#' @export
#' @rdname accessors
seasonal <- function(object) UseMethod("seasonal")

#' @export
#' @rdname accessors
trend <- function(object) UseMethod("trend")

#' @export
#' @rdname accessors
remainder <- function(object) UseMethod("remainder")

#' @export
#' @rdname accessors
getraw <- function(object) {
  object$data$raw
}

#' @export
#' @rdname accessors
remainder.stlplus <- function(object) {
  if(object$pars$fc.number == 0) {
    object$data$remainder
  } else {
    object$fc$remainder
  }
}

#' @export
#' @rdname accessors
fitted.stlplus <- function(object, ...) {
  if(object$pars$fc.number == 0) {
    object$data$seasonal + object$data$trend
  } else {
    object$data$seasonal + apply(object$fc[,1:object$pars$fc.number], 1, sum)
  }
}

#' @export
#' @rdname accessors
predict.stlplus <- function(object, ...) {
  if(object$pars$fc.number == 0) {
    object$data$seasonal + object$data$trend
  } else {
    object$data$seasonal + apply(object$fc[,1:object$pars$fc.number], 1, sum)
  }
}

#' @export
#' @rdname accessors
seasonal.stlplus <- function(object) {
  object$data$seasonal
}

#' @export
#' @rdname accessors
trend.stlplus <- function(object) {
  object$data$trend
}

# get frequency components:
#' @export
#' @rdname accessors
fc <- function(object, fcnum = 1) {
  if(! "stlplus" %in% class(object))
    stop("object not of class stlplus")

  if(is.null(object$fc))
    stop("there are no post-trend frequency components")

  if(fcnum > ncol(object$fc))
    stop("there are not that many frequency components")

  object$fc[,fcnum]
}

# setGeneric("time")
#' @export
#' @rdname accessors
time.stlplus <- function(x, ...) {
  if(length(x$t) == x$n) {
    x$t
  } else {
    c(1:x$n)
  }
}

# now for objects of class stl

#' @export
#' @rdname accessors
remainder.stl <- function(object) {
  as.numeric(object$time.series[,3])
}

#' @export
#' @rdname accessors
seasonal.stl <- function(object) {
  as.numeric(object$time.series[,1])
}

#' @export
#' @rdname accessors
trend.stl <- function(object) {
  as.numeric(object$time.series[,2])
}

# setGeneric("time")
#' @export
#' @rdname accessors
time.stl <- function(x, ...) {
  as.numeric(stats::time(x$time.series))
}

#' @export
#' @rdname accessors
predict.stl <- function(object, ...) {
  seasonal(object) + trend(object)
}

#' @export
#' @rdname accessors
fitted.stl <- function(object, ...) {
  seasonal(object) + trend(object)
}
