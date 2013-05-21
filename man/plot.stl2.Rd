\name{plot.stl2}
\alias{plot.stl2}
\title{Lattice plot of the raw, seasonal, trend, and remainder components}
\description{
Lattice plot of the raw, seasonal, trend, and remainder components.  If post-trend smoothing was done, these components will be plotted instead of the trend component.
}
\usage{
\method{plot}{stl2}(x, scales = list(y = list(relation = "sliced")), 
   type = "l", as.table = TRUE, strip = FALSE, strip.left = TRUE, 
   between = list(y = 0.5), layout = NULL, ...)
}
\arguments{
  \item{x}{object of class \code{"stl2"}.}
  \item{scales, type, as.table, strip, 
strip.left, between, layout, \dots}{parameters to be passed to xyplot.}
}
% \details{
% 
% }
\value{
object of class \code{"trellis"}.
}
\references{
R. B. Cleveland, W. S. Cleveland, J. E.  McRae, and I. Terpenning (1990) STL:  A  Seasonal-Trend  Decomposition  Procedure Based on Loess. \emph{Journal of Official Statistics}, \bold{6}, 3--73.
}
\author{Ryan Hafen}
% \note{ ~~further notes~~ 
%  ~Make other sections like Warning with \section{Warning }{\dots.} ~
% }
\seealso{\code{\link{stl2}}}
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--   or do  help(data=index)  for the standard data sets.
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
