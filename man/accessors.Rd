\name{accessors}
\alias{getraw}
\alias{seasonal}
\alias{trend}
\alias{remainder}
\alias{seasonal.stl}
\alias{seasonal.stl2}
\alias{trend.stl}
\alias{trend.stl2}
\alias{remainder.stl}
\alias{remainder.stl2}
\alias{time.stl}
\alias{time.stl2}
\alias{fitted.stl}
\alias{fitted.stl2}
\alias{predict.stl}
\alias{predict.stl2}
\alias{fc}

\title{Accessor functions for elements of an stl and stl2 object}
\description{
Retrieves the raw, seasonal, trend, remainder, or time components from an stl2 object.  The methods \code{seasonal.stl}, ... also exist as a convenience for extracting components from R's \code{stl()}.
}
\usage{
getraw(object)
\method{seasonal}{stl2}(object)
\method{trend}{stl2}(object)
\method{remainder}{stl2}(object)
\method{fitted}{stl2}(object, ...)
\method{time}{stl2}(x, ...)
fc(object, fcnum=1)

\method{seasonal}{stl}(object)
\method{trend}{stl}(object)
\method{remainder}{stl}(object)
\method{fitted}{stl}(object, ...)
\method{time}{stl}(x, ...)
}
\arguments{
  \item{x, object}{object of class \code{"stl"} or \code{"stl2"}.}
  \item{fcnum}{number of post-trend smoothing frequency component.}
  \item{...}{}
}
% \details{
% }
\value{
Returns a vector of either the \code{getraw} time series, the \code{seasonal}, \code{trend}, or \code{remainder} components, or the \code{time} values of the time series.  If \code{time}s are requested but were not supplied in the initial \code{stl2} call, the \code{1:n} vector is returned, where \code{n} is the number of data points.  The \code{fitted} method returns the sum of the seasonal and trend.
}
\references{
R. B. Cleveland, W. S. Cleveland, J.E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
}
\author{Ryan Hafen}
\note{
The \code{fitted} and \code{predict} methods are equivalent.  For objects of class \code{"stl2"}, these functions return the sum of all components but the remainder, including post-trend smoothing components.  Note also that the \code{trend} method for objects of class \code{"stl2"} only returns the trend component from the STL iterations, even when post-trend smoothing is done.
}
\seealso{\code{\link{stl2}}}
\examples{
co2.stl <- stl2(co2, t=as.vector(time(co2)), n.p=12, l.window=13, 
t.window=19, s.window=35, s.degree=1, sub.labels=substr(month.name, 1, 3))

plot(seasonal(co2.stl))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
